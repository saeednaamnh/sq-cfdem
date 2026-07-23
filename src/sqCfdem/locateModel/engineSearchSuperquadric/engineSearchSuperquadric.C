/* sqCfdem - engineSearchSuperquadric.C  (GPLv3) */

#include "error.H"
#include "engineSearchSuperquadric.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(engineSearchSuperquadric, 0);

addToRunTimeSelectionTable
(
    locateModel,
    engineSearchSuperquadric,
    dictionary
);

engineSearchSuperquadric::engineSearchSuperquadric
(
    const dictionary& dict,
    cfdemCloud& sm,
    word name
)
:
    engineSearch(dict, sm, "engine")   // reuse engineProps from the dict
{}

engineSearchSuperquadric::~engineSearchSuperquadric() {}

} // namespace Foam
