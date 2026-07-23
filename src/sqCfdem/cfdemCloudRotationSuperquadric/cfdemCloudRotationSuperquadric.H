/* sqCfdem - cfdemCloudRotationSuperquadric.H  (GPLv3)  [round 2]
   Concrete superquadric cloud. Name fixed by upstream declareCFDEMcloud.H/
   constructCFDEMcloud.H, which instantiate exactly this class when compiled
   with -DMS -Dsuperquadrics_flag.
   Overrides the base-class virtual accessors (which FatalError-abort for
   spheres): shape(), quat(), blockiness(), volume(), area(), shapeTypeName().
   Adds: quatComponents(), omega(), addHydroTorque() used by sqCfdem models. */
#ifndef cfdemCloudRotationSuperquadric_H
#define cfdemCloudRotationSuperquadric_H

#include "cfdemCloud.H"
#include "quaternion.H"

namespace Foam
{

class cfdemCloudRotationSuperquadric
:
    public cfdemCloud
{
    mutable int idShape_;
    mutable int idBlockinessLoc_;
    mutable int idQuatLoc_;
    mutable int idOmega_;
    mutable int idVolumeLoc_;
    mutable int idAreaLoc_;
    mutable int idHdTorque_;

public:
    cfdemCloudRotationSuperquadric(const fvMesh& mesh);
    ~cfdemCloudRotationSuperquadric();

    // --- data exchange ------------------------------------------------
    virtual void getDEMdata();
    virtual void giveDEMdata();

    // --- overrides of base virtuals (base aborts) ---------------------
    virtual vector     shape(int i) const;
    virtual vector2D   blockiness(int i) const;
    virtual quaternion quat(int i) const;
    virtual scalar     volume(int i) const;
    virtual scalar     area(int i) const;
    virtual inline word shapeTypeName() { return "superquadric"; }

    // --- convenience for sqCfdem models -------------------------------
    void quatComponents(int i, scalar& w, scalar& x,
                        scalar& y, scalar& z) const;
    vector omegaSQ(int i) const;
    void addHydroTorque(int i, const vector& T) const;
};

} // namespace Foam
#endif
