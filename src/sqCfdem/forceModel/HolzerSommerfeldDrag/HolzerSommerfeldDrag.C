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
    useDiFeliceCorrection_
    (
        propsDict_.lookupOrDefault<Switch>("voidageCorrection", true)
    )
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

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        const label cellI =
            static_cast<label>(particleCloud_.cfdemCloud::cellIDs()[index][0]);
        if (cellI < 0) continue;

        const vector Ufluid = U_[cellI];
        const vector Us = particleCloud_.velocity(index);
        const vector Ur = Ufluid - Us;
        const scalar magUr = mag(Ur);
        if (magUr < SMALL) continue;

        const vector shp = particleCloud_.shape(index);
        const scalar a = shp.x();
        const scalar c = shp.z();

        const quaternion q = particleCloud_.quat(index);
        vector k = sqGeom::symmetryAxis(q.w(), q.v().x(), q.v().y(), q.v().z());
        k /= mag(k) + SMALL;

        const scalar cosA = mag(k & (Ur/magUr));

        const scalar rho  = rhoField[cellI];
        const scalar nuf  = nufField[cellI];
        const scalar dpEq = sqGeom::dEq(a,c);
        const scalar Re   = dpEq*magUr/nuf;
        const scalar Cd   = sqGeom::CdHolzerSommerfeld(Re, a, c, cosA);

        scalar fEps = 1.0;
        if (useDiFeliceCorrection_)
        {
            const scalar eps = voidfraction_[cellI];
            fEps = pow(eps, -sqGeom::diFeliceBeta(Re));
        }

        const scalar Aperp = sqGeom::projectedArea(a, c, cosA);
        const scalar dragCoeff = 0.5*rho*Cd*Aperp*magUr*fEps;  // F = coeff*Ur
        const vector drag = dragCoeff*Ur;

        forceSubM(0).partToArray(index, drag, vector::zero, Ufluid, dragCoeff);
    }
}

} // namespace Foam
