/* sqCfdem - engineSearchSuperquadric.H  (GPLv3)
   Locate model for superquadric particles. Particle *centres* are located
   exactly as for spheres (shape extent is handled by the satellite points
   in dividedVoidFractionSuperquadric), so this class simply reuses
   engineSearch under the suffixed name that upstream's newLocateModel
   dispatch requests when particleShapeType is superquadric. */
#ifndef engineSearchSuperquadric_H
#define engineSearchSuperquadric_H

#include "engineSearch.H"

namespace Foam
{

class engineSearchSuperquadric
:
    public engineSearch
{
public:
    TypeName("engineSuperquadric");

    engineSearchSuperquadric
    (
        const dictionary& dict,
        cfdemCloud& sm,
        word name = word("")
    );

    ~engineSearchSuperquadric();
};

} // namespace Foam
#endif
