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

    // unit-ball quasi-uniform marker set (Fibonacci shells + centre),
    // scaled per particle by semi-axes (a,b,c) at use time.
    unitMarkers_.setSize(nSat_);
    unitMarkers_[0] = vector::zero;
    const label nShell = nSat_ - 1;
    const scalar golden = M_PI*(3.0 - std::sqrt(5.0));
    for (label i = 0; i < nShell; i++)
    {
        const scalar r  = std::cbrt((i + 0.5)/nShell);
        const scalar zz = 1.0 - 2.0*(i + 0.5)/nShell;
        const scalar rho = std::sqrt(std::max(0.0, 1.0 - zz*zz));
        const scalar th = golden*i;
        unitMarkers_[i+1] = r*vector(rho*std::cos(th), rho*std::sin(th), zz);
    }
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
        const scalar Vp    = particleCloud_.volume(index);
        const scalar subV  = Vp/nSat_;

        const quaternion q = particleCloud_.quat(index);
        const scalar qw = q.w();
        const vector qv = q.v();
        const tensor R  = sqGeom::quatToRot(qw, qv.x(), qv.y(), qv.z());

        label nSlots = 0;
        labelList slotCell(maxCellsPerParticle_, label(-1));
        label lostSat = 0;

        for (label s = 0; s < nSat_; s++)
        {
            vector xb
            (
                shp.x()*unitMarkers_[s].x(),
                shp.y()*unitMarkers_[s].y(),
                shp.z()*unitMarkers_[s].z()
            );
            vector xw  = pos + (R & xb);   // non-const: locate API takes refs
            label seed = cellCentre;

            const label cellI =
                particleCloud_.locateM().findSingleCell(xw, seed);

            if (cellI < 0) { lostSat++; continue; }  // TODO(build-2): stash +
                                                     // parallel scatter-gather

            label slot = -1;
            for (label t = 0; t < nSlots; t++)
                if (slotCell[t] == cellI) { slot = t; break; }
            if (slot < 0)
            {
                if (nSlots >= maxCellsPerParticle_) slot = 0;
                else
                {
                    slot = nSlots++;
                    slotCell[slot] = cellI;
                    particleCloud_.cfdemCloud::cellIDs()[index][slot] = cellI;
                }
            }

            particleWeights[index][slot] += 1.0/nSat_;
            particleVolumes[index][slot] += subV;
            voidfractionNext_[cellI]     -= subV/cellVol[cellI];
        }
        cellsPerParticle_[index][0] = scalar(max(nSlots, label(1)));

        if (lostSat > 0 && lostSat < nSat_)
        {
            const scalar f = scalar(nSat_)/scalar(nSat_ - lostSat);
            for (label t = 0; t < nSlots; t++)
            {
                particleWeights[index][t] *= f;
                particleVolumes[index][t] *= f;
            }
        }

        particleV[index][0] = Vp;
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
