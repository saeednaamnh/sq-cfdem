/*---------------------------------------------------------------------------*\
    sqCfdem - dividedVoidFractionSuperquadric   (GPLv3)  [build-fix round 2]
    Satellite-point void fraction for prolate spheroids.
    Registered as "dividedSuperquadric" (upstream <model>+"Superquadric"
    dispatch when particleShapeType is superquadric).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "dividedVoidFractionSuperquadric.H"
#include "addToRunTimeSelectionTable.H"
#include "locateModel.H"
#include "quaternion.H"
#include "spheroidGeometry.H"

namespace Foam
{

defineTypeNameAndDebug(dividedVoidFractionSuperquadric, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    dividedVoidFractionSuperquadric,
    dictionary
);

dividedVoidFractionSuperquadric::dividedVoidFractionSuperquadric
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict("dividedSuperquadricProps")),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
    nSat_(propsDict_.lookupOrDefault<label>("nSatellites", 64)),
    weight_(propsDict_.lookupOrDefault<scalar>("weight", 1.0))
{
    maxCellsPerParticle_ =
        propsDict_.lookupOrDefault<label>("maxCellsPerParticle", 30);

    // Unit-ball markers in antipodal pairs (centroid exactly zero), scaled per
    // particle by the semi-axes (a,b,c) at use time; see sqGeom::unitBallMarkers.
    // (The previous set tied radius to direction, which shifted each particle's
    // deposited volume by about 11 % of c along its own axis.)
    const std::vector<vector> m = sqGeom::unitBallMarkers(nSat_);
    unitMarkers_.setSize(label(m.size()));
    forAll(unitMarkers_, i) unitMarkers_[i] = m[i];
    Info << type() << ": " << nSat_ << " satellite points per particle" << endl;
}

dividedVoidFractionSuperquadric::~dividedVoidFractionSuperquadric() {}

void dividedVoidFractionSuperquadric::setvoidFraction
(
    double** const& mask,
    double**& voidfractions,
    double**& particleWeights,
    double**& particleVolumes,
    double**& particleV
) const
{
    reAllocArrays();

    voidfractionNext_.primitiveFieldRef() = 1.0;

    const scalarField& cellVol = voidfractionNext_.mesh().V();
    static label nCalls = 0;  ++nCalls;
    label lostThisCall = 0, overflowThisCall = 0;

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        for (int s = 0; s < maxCellsPerParticle_; s++)
        {
            particleWeights[index][s] = 0.;
            particleVolumes[index][s] = 0.;
        }
        cellsPerParticle_[index][0] = 1.;

        const label cellCentre =
            static_cast<label>(particleCloud_.cfdemCloud::cellIDs()[index][0]);
        if (cellCentre < 0) continue;

        const vector pos   = particleCloud_.position(index);
        const vector shp   = particleCloud_.shape(index);      // (a,b,c)
        const vector2D blk = particleCloud_.blockiness(index);
        const std::string why =
            sqGeom::unsupportedShape(shp.x(), shp.y(), shp.z(), blk.x(), blk.y());
        if (!why.empty())
        {
            FatalError << "dividedSuperquadric: particle " << index << ": "
                       << why.c_str() << abort(FatalError);
        }
        const scalar Vp    = particleCloud_.volume(index);
        const scalar subV  = Vp/unitMarkers_.size();

        const quaternion q = particleCloud_.quat(index);
        const scalar qw = q.w();
        const vector qv = q.v();
        const tensor R  = sqGeom::quatToRot(qw, qv.x(), qv.y(), qv.z());

        // Slot 0 is always the cell of the particle centre: the drag and torque
        // models read the fluid there (cellIDs[index][0]).
        label nSlots = 1;
        labelList slotCell(maxCellsPerParticle_, label(-1));
        slotCell[0] = cellCentre;

        for (label s = 0; s < label(unitMarkers_.size()); s++)
        {
            vector xb
            (
                shp.x()*unitMarkers_[s].x(),
                shp.y()*unitMarkers_[s].y(),
                shp.z()*unitMarkers_[s].z()
            );
            vector xw  = pos + (R & xb);   // non-const: locate API takes refs
            label seed = cellCentre;

            label cellI = particleCloud_.locateM().findSingleCell(xw, seed);

            // A marker outside this processor's cells (another processor, or
            // outside the mesh) is deposited in the centre cell, so that no
            // particle volume is lost; its location is then approximate.
            if (cellI < 0) { cellI = cellCentre; lostThisCall++; }

            label slot = -1;
            for (label t = 0; t < nSlots; t++)
                if (slotCell[t] == cellI) { slot = t; break; }
            if (slot < 0)
            {
                if (nSlots >= maxCellsPerParticle_) { slot = 0; overflowThisCall++; }
                else
                {
                    slot = nSlots++;
                    slotCell[slot] = cellI;
                    particleCloud_.cfdemCloud::cellIDs()[index][slot] = cellI;
                }
            }

            const label cellS = slotCell[slot];
            particleWeights[index][slot] += 1.0/unitMarkers_.size();
            particleVolumes[index][slot] += subV;
            voidfractionNext_[cellS]     -= subV/cellVol[cellS];
        }
        cellsPerParticle_[index][0] = scalar(nSlots);

        particleV[index][0] = Vp;
    }

    if ((lostThisCall > 0 || overflowThisCall > 0) && (nCalls <= 3 || nCalls % 1000 == 0))
    {
        Pout<< type() << ": call " << nCalls << ": " << lostThisCall
            << " markers outside this processor's cells and " << overflowThisCall
            << " beyond maxCellsPerParticle were deposited in the particle's"
            << " centre cell (volume conserved, location approximate)" << endl;
    }

    forAll(voidfractionNext_, cI)
        voidfractionNext_[cI] = max(voidfractionNext_[cI], alphaMin_);

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        const int nC = static_cast<int>(cellsPerParticle_[index][0]);
        for (int s = 0; s < nC; s++)
        {
            const label cellI =
                static_cast<label>
                (particleCloud_.cfdemCloud::cellIDs()[index][s]);
            if (cellI >= 0)
                voidfractions[index][s] = voidfractionNext_[cellI];
        }
    }
}

} // namespace Foam
