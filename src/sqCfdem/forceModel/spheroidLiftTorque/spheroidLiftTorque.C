/* sqCfdem - spheroidLiftTorque.C  (GPLv3)  [round 2]
   Lift + pitching torque. Torque push-back goes through the derived
   cloud's addHydroTorque(); a dynamic_cast at construction guarantees
   the solver actually instantiated the superquadric cloud. */

#include "error.H"
#include "spheroidLiftTorque.H"
#include "addToRunTimeSelectionTable.H"
#include "quaternion.H"
#include "spheroidGeometry.H"
#include "cfdemCloudRotationSuperquadric.H"

namespace Foam
{

defineTypeNameAndDebug(spheroidLiftTorque, 0);
addToRunTimeSelectionTable(forceModel, spheroidLiftTorque, dictionary);

spheroidLiftTorque::spheroidLiftTorque(const dictionary& dict, cfdemCloud& sm)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField>(velFieldName_)),
    f1_(readScalar(propsDict_.lookup("liftPrefactor"))),
    f2_(readScalar(propsDict_.lookup("liftReExponent"))),
    g1_(readScalar(propsDict_.lookup("torquePrefactor"))),
    g2_(readScalar(propsDict_.lookup("torqueReExponent"))),
    spinDamping_(propsDict_.lookupOrDefault<Switch>("spinDamping", true)),
    sqCloud_(dynamic_cast<cfdemCloudRotationSuperquadric*>(&sm))
{
    if (!sqCloud_)
    {
        FatalError << "spheroidLiftTorque requires the superquadric cloud "
                   << "(use cfdemSolverPisoSQ, not a sphere solver)."
                   << abort(FatalError);
    }
    setForceSubModels(propsDict_);
    forceSubM(0).readSwitches();
}

spheroidLiftTorque::~spheroidLiftTorque() {}

void spheroidLiftTorque::setForce() const
{
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        const label cellI =
            static_cast<label>(particleCloud_.cfdemCloud::cellIDs()[index][0]);
        if (cellI < 0) continue;

        const vector Ufluid = U_[cellI];
        const vector Ur = Ufluid - particleCloud_.velocity(index);
        const scalar magUr = mag(Ur);
        if (magUr < SMALL) continue;
        const vector uhat = Ur/magUr;

        const vector shp = particleCloud_.shape(index);
        const scalar a = shp.x(), c = shp.z();

        const quaternion q = particleCloud_.quat(index);
        vector k = sqGeom::symmetryAxis(q.w(), q.v().x(), q.v().y(), q.v().z());
        k /= mag(k)+SMALL;

        const scalar cosA = k & uhat;
        const vector kxu  = k ^ uhat;
        const scalar sinA = mag(kxu);
        if (sinA < SMALL) continue;      // axis-aligned: no lift/torque

        const scalar rho  = rhoField[cellI];
        const scalar nuf  = nufField[cellI];
        const scalar dpEq = sqGeom::dEq(a,c);
        const scalar Re   = dpEq*magUr/nuf;
        const scalar Aref = 0.25*M_PI*dpEq*dpEq;

        vector eL = (kxu ^ uhat); eL /= mag(eL)+SMALL;
        const scalar Cl =
            sqGeom::ClSpheroid(Re, sinA, mag(cosA), f1_, f2_);
        const vector lift = 0.5*rho*Cl*Aref*magUr*magUr*eL*sign(cosA);

        const vector eT = kxu/(mag(kxu)+SMALL);
        const scalar Ct =
            sqGeom::CtSpheroid(Re, sinA, mag(cosA), g1_, g2_);
        const vector Tp =
            0.5*rho*Ct*Aref*(0.5*dpEq)*magUr*magUr*eT*sign(cosA);

        forceSubM(0).partToArray(index, lift, vector::zero, Ufluid, scalar(0));
        sqCloud_->addHydroTorque(index, Tp);
    }
}

} // namespace Foam
