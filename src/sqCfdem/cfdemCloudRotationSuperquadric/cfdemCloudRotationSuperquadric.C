/* sqCfdem - cfdemCloudRotationSuperquadric.C  (GPLv3)  [round 2] */

#include "cfdemCloudRotationSuperquadric.H"

namespace Foam
{

cfdemCloudRotationSuperquadric::cfdemCloudRotationSuperquadric
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    idShape_(-1), idBlockinessLoc_(-1), idQuatLoc_(-1),
    idOmega_(-1), idVolumeLoc_(-1), idAreaLoc_(-1), idHdTorque_(-1)
{
    // pull=true : DEM -> CFD  (names/types mirror fix_cfd_coupling_force.cpp)
    registerFieldsToDEM("shape",      "vector-atom",     idShape_,        true);
    registerFieldsToDEM("blockiness", "vector2D-atom",   idBlockinessLoc_,true);
    registerFieldsToDEM("quaternion", "quaternion-atom", idQuatLoc_,      true);
    registerFieldsToDEM("omega",      "vector-atom",     idOmega_,        true);
    registerFieldsToDEM("volume",     "scalar-atom",     idVolumeLoc_,    true);
    registerFieldsToDEM("area",       "scalar-atom",     idAreaLoc_,      true);
    // push : CFD -> DEM (hydrodynamic torque; LIGGGHTS adds to torque[i])
    registerFieldsToDEM("hdtorque",   "vector-atom",     idHdTorque_,     false);
}

cfdemCloudRotationSuperquadric::~cfdemCloudRotationSuperquadric() {}

void cfdemCloudRotationSuperquadric::getDEMdata()
{
    cfdemCloud::getDEMdata();
}

void cfdemCloudRotationSuperquadric::giveDEMdata()
{
    cfdemCloud::giveDEMdata();
}

vector cfdemCloudRotationSuperquadric::shape(int i) const
{
    return vector(fieldsToDEM[idShape_][i][0],
                  fieldsToDEM[idShape_][i][1],
                  fieldsToDEM[idShape_][i][2]);
}

vector2D cfdemCloudRotationSuperquadric::blockiness(int i) const
{
    return vector2D(fieldsToDEM[idBlockinessLoc_][i][0],
                    fieldsToDEM[idBlockinessLoc_][i][1]);
}

quaternion cfdemCloudRotationSuperquadric::quat(int i) const
{
    return quaternion
    (
        fieldsToDEM[idQuatLoc_][i][0],
        vector
        (
            fieldsToDEM[idQuatLoc_][i][1],
            fieldsToDEM[idQuatLoc_][i][2],
            fieldsToDEM[idQuatLoc_][i][3]
        )
    );
}

scalar cfdemCloudRotationSuperquadric::volume(int i) const
{
    return fieldsToDEM[idVolumeLoc_][i][0];
}

scalar cfdemCloudRotationSuperquadric::area(int i) const
{
    return fieldsToDEM[idAreaLoc_][i][0];
}

void cfdemCloudRotationSuperquadric::quatComponents
(
    int i, scalar& w, scalar& x, scalar& y, scalar& z
) const
{
    w = fieldsToDEM[idQuatLoc_][i][0];
    x = fieldsToDEM[idQuatLoc_][i][1];
    y = fieldsToDEM[idQuatLoc_][i][2];
    z = fieldsToDEM[idQuatLoc_][i][3];
}

vector cfdemCloudRotationSuperquadric::omegaSQ(int i) const
{
    return vector(fieldsToDEM[idOmega_][i][0],
                  fieldsToDEM[idOmega_][i][1],
                  fieldsToDEM[idOmega_][i][2]);
}

void cfdemCloudRotationSuperquadric::addHydroTorque
(
    int i, const vector& T
) const
{
    fieldsToDEM[idHdTorque_][i][0] += T.x();
    fieldsToDEM[idHdTorque_][i][1] += T.y();
    fieldsToDEM[idHdTorque_][i][2] += T.z();
}

} // namespace Foam
