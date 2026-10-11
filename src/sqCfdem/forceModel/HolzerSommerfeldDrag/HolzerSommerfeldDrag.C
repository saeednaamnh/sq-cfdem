/* sqCfdem - HolzerSommerfeldDrag.C  (GPLv3)  [round 2]
   Orientation-dependent drag for prolate spheroids.
   Uses base-class virtuals shape()/quat() (implemented by
   cfdemCloudRotationSuperquadric). Cell-centre fluid velocity for now;
   cellPoint interpolation is TODO(build-3). */

#include "error.H"
#include "HolzerSommerfeldDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "quaternion.H"
#include "spheroidGeometry.H"
#include "cfdemCloudRotationSuperquadric.H"

namespace Foam
{

defineTypeNameAndDebug(HolzerSommerfeldDrag, 0);
addToRunTimeSelectionTable(forceModel, HolzerSommerfeldDrag, dictionary);

HolzerSommerfeldDrag::HolzerSommerfeldDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField>(velFieldName_)),
    voidfractionFieldName_
    (
        propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")
    ),
    voidfraction_
    (
        sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)
    ),
    // No default: whether the drag is corrected for the local solid fraction
    // changes the physics of every dense case, so the user must state it.
    useDiFeliceCorrection_(propsDict_.lookup("voidageCorrection"))
{
    setForceSubModels(propsDict_);
    forceSubM(0).readSwitches();
    particleCloud_.checkCG(true);
}

HolzerSommerfeldDrag::~HolzerSommerfeldDrag() {}

void HolzerSommerfeldDrag::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    const cfdemCloudRotationSuperquadric* sq =
        dynamic_cast<const cfdemCloudRotationSuperquadric*>(&particleCloud_);
    if (!sq)
        FatalError << "HolzerSommerfeldDrag requires the superquadric cloud"
                   << " (cfdemSolverPisoSQ)." << abort(FatalError);

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        if (sq->excluded(index)) continue;      // e.g. wall spheres: no drag exchange
        const label cellI =
            static_cast<label>(particleCloud_.cfdemCloud::cellIDs()[index][0]);
        if (cellI < 0) continue;

        const vector Ufluid = U_[cellI];
        const vector Us = particleCloud_.velocity(index);
        const vector Ur = Ufluid - Us;
        const scalar magUr = mag(Ur);
        if (magUr < SMALL) continue;

        const vector shp = particleCloud_.shape(index);
        const vector2D blk = particleCloud_.blockiness(index);
        const std::string why =
            sqGeom::unsupportedShape(shp.x(), shp.y(), shp.z(), blk.x(), blk.y());
        if (!why.empty())
        {
            FatalError << "HolzerSommerfeldDrag: particle " << index << ": "
                       << why.c_str() << abort(FatalError);
        }
        const scalar a = shp.x();
        const scalar c = shp.z();

        const quaternion q = particleCloud_.quat(index);
        vector k = sqGeom::symmetryAxis(q.w(), q.v().x(), q.v().y(), q.v().z());
        k /= mag(k) + SMALL;

        const scalar cosA = mag(k & (Ur/magUr));

        // Hoelzer & Sommerfeld C_D referred to the volume-equivalent sphere,
        // optionally with the Di Felice voidage correction in its standard
        // (superficial-velocity) form; see sqGeom::dragCoefficient.
        scalar Re = 0;
        const scalar dragCoeff = sqGeom::dragCoefficient
        (
            rhoField[cellI], nufField[cellI], a, c, cosA, magUr,
            voidfraction_[cellI], useDiFeliceCorrection_, Re
        );                                                   // F = coeff*Ur
        const vector drag = dragCoeff*Ur;

        forceSubM(0).partToArray(index, drag, vector::zero, Ufluid, dragCoeff);
    }
}

} // namespace Foam
