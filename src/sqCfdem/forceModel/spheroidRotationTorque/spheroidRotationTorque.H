/* sqCfdem - spheroidRotationTorque.H  (GPLv3)
   Exact-Stokes rotational torque on a prolate spheroid:
     T = 8 pi mu a^2 c [ XC dW_par + YC dW_perp + YH k x (E.k) ],
   dW = 0.5 curl(U) - omega_p, E = symm(grad U). Drives Jeffery dynamics
   (V3). Torque-only model: contributes no force to Ksl.               */
#ifndef spheroidRotationTorque_H
#define spheroidRotationTorque_H

#include "forceModel.H"

namespace Foam
{

class cfdemCloudRotationSuperquadric;

class spheroidRotationTorque
:
    public forceModel
{
    dictionary            propsDict_;
    word                  velFieldName_;
    const volVectorField& U_;
    mutable cfdemCloudRotationSuperquadric* sqCloud_;

public:
    TypeName("spheroidRotationTorque");

    spheroidRotationTorque(const dictionary& dict, cfdemCloud& sm);
    ~spheroidRotationTorque();

    void setForce() const;
};

} // namespace Foam
#endif
