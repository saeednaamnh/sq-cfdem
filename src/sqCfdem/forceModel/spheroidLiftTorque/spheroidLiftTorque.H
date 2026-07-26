/* sqCfdem - spheroidLiftTorque.H  (GPLv3)  [round 2] */
#ifndef spheroidLiftTorque_H
#define spheroidLiftTorque_H

#include "forceModel.H"

namespace Foam
{

class cfdemCloudRotationSuperquadric;

class spheroidLiftTorque
:
    public forceModel
{
    dictionary            propsDict_;
    word                  velFieldName_;
    const volVectorField& U_;
    scalar f1_, f2_;   // lift prefactor/exponent
    scalar g1_, g2_;   // pitching torque prefactor/exponent
    Switch spinDamping_;
    mutable cfdemCloudRotationSuperquadric* sqCloud_;

public:
    TypeName("spheroidLiftTorque");

    spheroidLiftTorque(const dictionary& dict, cfdemCloud& sm);
    ~spheroidLiftTorque();

    void setForce() const;
};

} // namespace Foam
#endif
