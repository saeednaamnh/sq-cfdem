/* sqCfdem - dividedVoidFractionSuperquadric.H  (GPLv3) */
#ifndef dividedVoidFractionSuperquadric_H
#define dividedVoidFractionSuperquadric_H

#include "voidFractionModel.H"

namespace Foam
{

class dividedVoidFractionSuperquadric
:
    public voidFractionModel
{
    dictionary   propsDict_;
    scalar       alphaMin_;
    label        nSat_;
    scalar       weight_;
    List<vector> unitMarkers_;   // body-frame unit-ball satellite lattice

public:
    // NOTE: name feeds upstream dispatch <model>+"Superquadric"
    TypeName("dividedSuperquadric");

    dividedVoidFractionSuperquadric(const dictionary& dict, cfdemCloud& sm);
    ~dividedVoidFractionSuperquadric();

    void setvoidFraction
    (
        double** const& mask,
        double**& voidfractions,
        double**& particleWeights,
        double**& particleVolumes,
        double**& particleV
    ) const;
};

} // namespace Foam
#endif
