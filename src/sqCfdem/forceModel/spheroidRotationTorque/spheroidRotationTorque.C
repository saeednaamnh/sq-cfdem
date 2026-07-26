/* sqCfdem - spheroidRotationTorque.C  (GPLv3) */

#include "error.H"
#include "spheroidRotationTorque.H"
#include "addToRunTimeSelectionTable.H"
#include "quaternion.H"
#include "fvCFD.H"
#include "spheroidGeometry.H"
#include "cfdemCloudRotationSuperquadric.H"

namespace Foam
{

defineTypeNameAndDebug(spheroidRotationTorque, 0);
addToRunTimeSelectionTable(forceModel, spheroidRotationTorque, dictionary);

spheroidRotationTorque::spheroidRotationTorque
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField>(velFieldName_)),
    sqCloud_(nullptr)   // cast deferred: cloud not fully constructed yet
{
    setForceSubModels(propsDict_);
    forceSubM(0).readSwitches();
}

spheroidRotationTorque::~spheroidRotationTorque() {}

void spheroidRotationTorque::setForce() const
{
    if (!sqCloud_)
    {
        sqCloud_ = dynamic_cast<cfdemCloudRotationSuperquadric*>(&particleCloud_);
        if (!sqCloud_)
            FatalError << "spheroidRotationTorque requires the superquadric "
                       << "cloud (cfdemSolverPisoSQ)." << abort(FatalError);
    }
    const volScalarField& nufField = forceSubM(0).nuField();
    const volScalarField& rhoField = forceSubM(0).rhoField();

    // velocity-gradient fields once per call
    const volTensorField gradU(fvc::grad(U_));
    const volVectorField curlU(fvc::curl(U_));

    for (int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
        const label cellI =
            static_cast<label>(particleCloud_.cfdemCloud::cellIDs()[index][0]);
        if (cellI < 0) continue;

        const vector shp = particleCloud_.shape(index);
        const scalar a = shp.x(), c = shp.z();

        const quaternion q = particleCloud_.quat(index);
        vector k = sqGeom::symmetryAxis(q.w(), q.v().x(), q.v().y(), q.v().z());
        k /= mag(k) + SMALL;

        const vector Wf = 0.5*curlU[cellI];
        const vector wp = sqCloud_->omegaSQ(index);
        const vector dW = Wf - wp;

        const symmTensor E = symm(gradU[cellI]);
        const vector Ek = E & k;

        scalar XC, YC, YH;
        sqGeom::rotResist(a, c, XC, YC, YH);

        const scalar mu = rhoField[cellI]*nufField[cellI];
        const scalar pref = 8.0*M_PI*mu*a*a*c;

        const vector dWpar  = (k & dW)*k;
        const vector dWperp = dW - dWpar;

        const vector T = pref*( XC*dWpar + YC*dWperp + YH*(k ^ Ek) );

        sqCloud_->addHydroTorque(index, T);
    }
}

} // namespace Foam
