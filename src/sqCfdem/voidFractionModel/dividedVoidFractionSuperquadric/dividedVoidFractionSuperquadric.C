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
#include "cfdemCloudRotationSuperquadric.H"

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
    weight_(propsDict_.lookupOrDefault<scalar>("weight", 1.0)),
    // No default: which directions are periodic decides whether a marker beyond
    // the domain is wrapped to the opposite side or lies outside the fluid.
    periodic_(propsDict_.lookup("periodicDirections")),
    domain_(sm.mesh().points(), true)
{
    for (direction d = 0; d < vector::nComponents; d++)
    {
        if (periodic_[d] != 0 && periodic_[d] != 1)
        {
            FatalError << "dividedSuperquadric: periodicDirections must hold 0 or 1"
                       << " per direction, e.g. (1 1 0); found " << periodic_
                       << abort(FatalError);
        }
    }
    Info<< "dividedSuperquadric: fluid domain " << domain_ << ", periodic "
        << periodic_ << endl;
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
    label lostThisCall = 0, overflowThisCall = 0, wrappedThisCall = 0, outsideThisCall = 0;
    scalar Vtotal = 0, Vinside = 0;
    label nExcluded = 0;
    const cfdemCloudRotationSuperquadric* sq =
        dynamic_cast<const cfdemCloudRotationSuperquadric*>(&particleCloud_);
    if (!sq)
        FatalError << "dividedSuperquadric requires the superquadric cloud"
                   << " (cfdemSolverPisoSQ)." << abort(FatalError);
    const scalar lo[3] = {domain_.min().x(), domain_.min().y(), domain_.min().z()};
    const scalar hi[3] = {domain_.max().x(), domain_.max().y(), domain_.max().z()};
    const bool per[3] = {periodic_.x() > 0.5, periodic_.y() > 0.5, periodic_.z() > 0.5};

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        for (int s = 0; s < maxCellsPerParticle_; s++)
        {
            particleWeights[index][s] = 0.;
            particleVolumes[index][s] = 0.;
        }
        cellsPerParticle_[index][0] = 1.;
        particleV[index][0] = 0;


        // particles not coupled to the fluid (e.g. wall spheres): no volume, no forces.
        // Counted only where the particle centre is local (every processor holds all particles).
        if (sq->excluded(index))
        {
            if (particleCloud_.cfdemCloud::cellIDs()[index][0] >= 0) nExcluded++;
            continue;
        }

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

            // Periodic directions: wrap into the box. Wall directions: a marker
            // beyond the domain is outside the fluid, and so is its volume.
            scalar xp[3] = {xw.x(), xw.y(), xw.z()};
            const int place = sqGeom::placeMarker(xp, lo, hi, per);
            if (place == sqGeom::markerOutside) { outsideThisCall++; continue; }
            if (place == sqGeom::markerWrapped) { wrappedThisCall++; xw = vector(xp[0], xp[1], xp[2]); }

            label seed = cellCentre;
            label cellI = particleCloud_.locateM().findSingleCell(xw, seed);

            // Inside the fluid but in another processor's cells: deposited in the
            // centre cell, so that no volume is lost (location approximate).
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
            particleVolumes[index][slot] += subV;
            voidfractionNext_[cellS]     -= subV/cellVol[cellS];
        }
        cellsPerParticle_[index][0] = scalar(nSlots);

        // weights: share of the particle's volume inside the fluid, per cell
        scalar Vin = 0;
        for (label t = 0; t < nSlots; t++) Vin += particleVolumes[index][t];
        if (Vin > 0)
            for (label t = 0; t < nSlots; t++) particleWeights[index][t] = particleVolumes[index][t]/Vin;
        else
            particleWeights[index][0] = 1.0;
        Vtotal  += Vp;
        Vinside += Vin;

        particleV[index][0] = Vp;
    }

    // global volume budget (every processor takes part in the reductions)
    reduce(Vtotal, sumOp<scalar>());           reduce(Vinside, sumOp<scalar>());
    reduce(lostThisCall, sumOp<label>());      reduce(overflowThisCall, sumOp<label>());
    reduce(wrappedThisCall, sumOp<label>());   reduce(outsideThisCall, sumOp<label>());
    reduce(nExcluded, sumOp<label>());
    if (nCalls <= 3 || nCalls % 1000 == 0)
    {
        Info<< type() << ": call " << nCalls << ": " << nExcluded
            << " particles not coupled (couplingFlag 0); volume budget: particles "
            << Vtotal << ", inside the fluid domain " << Vinside << " ("
            << (Vtotal > 0 ? 100.0*(Vtotal - Vinside)/Vtotal : 0.0)
            << " % outside, e.g. wall spheres beyond the walls); markers wrapped"
            << " across periodic boundaries " << wrappedThisCall << ", outside the"
            << " domain " << outsideThisCall << ", on another processor (deposited in"
            << " the centre cell) " << lostThisCall << ", beyond maxCellsPerParticle "
            << overflowThisCall << endl;
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
