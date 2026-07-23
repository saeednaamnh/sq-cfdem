/* sqCfdem - HolzerSommerfeldDrag.H  (GPLv3) */
#ifndef HolzerSommerfeldDrag_H
#define HolzerSommerfeldDrag_H

#include "forceModel.H"

namespace Foam
{

class HolzerSommerfeldDrag
:
    public forceModel
{
    dictionary          propsDict_;
    word                velFieldName_;
    const volVectorField& U_;
    word                voidfractionFieldName_;
    const volScalarField& voidfraction_;
    Switch              useDiFeliceCorrection_;

public:
    TypeName("HolzerSommerfeldDrag");

    HolzerSommerfeldDrag(const dictionary& dict, cfdemCloud& sm);
    ~HolzerSommerfeldDrag();

    void setForce() const;
};

} // namespace Foam
#endif
