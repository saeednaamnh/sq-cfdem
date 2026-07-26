/*---------------------------------------------------------------------------*\
    sqCfdem - superquadric extension for CFDEMcoupling-PUBLIC
    License: GPLv3 (inherited from CFDEMcoupling)

    spheroidGeometry.H
    Header-only helpers for prolate spheroids (a = b < c), n1 = n2 = 2.
    All functions are pure math: no OpenFOAM state, unit-agnostic.
\*---------------------------------------------------------------------------*/

#ifndef spheroidGeometry_H
#define spheroidGeometry_H

#include "vector.H"
#include "tensor.H"
#include <cmath>

namespace Foam
{
namespace sqGeom
{

//- Quaternion (w,x,y,z) -> rotation matrix (body->world).
//  LIGGGHTS superquadric convention: symmetry axis = third column R.col(2).
inline tensor quatToRot(const scalar w, const scalar x,
                        const scalar y, const scalar z)
{
    return tensor
    (
        1-2*(y*y+z*z),  2*(x*y-w*z),    2*(x*z+w*y),
        2*(x*y+w*z),    1-2*(x*x+z*z),  2*(y*z-w*x),
        2*(x*z-w*y),    2*(y*z+w*x),    1-2*(x*x+y*y)
    );
}

//- Symmetry axis k (world frame) from quaternion
inline vector symmetryAxis(const scalar w, const scalar x,
                           const scalar y, const scalar z)
{
    // third column of R
    return vector(2*(x*z+w*y), 2*(y*z-w*x), 1-2*(x*x+y*y));
}

//- Volume-equivalent diameter, V = 4/3 pi a^2 c
inline scalar dEq(const scalar a, const scalar c)
{
    return 2.0*std::cbrt(a*a*c);
}

//- Surface area of a prolate spheroid (a=b<c), exact
inline scalar surfaceArea(const scalar a, const scalar c)
{
    if (c <= a*(1.0+1e-12))          // sphere fallback
        return 4.0*M_PI*a*a;
    const scalar e = std::sqrt(1.0 - (a*a)/(c*c));   // eccentricity
    return 2.0*M_PI*a*a*(1.0 + (c/(a*e))*std::asin(e));
}

//- Sphericity Phi = A_sphere(dEq)/A_particle
inline scalar sphericity(const scalar a, const scalar c)
{
    const scalar d = dEq(a,c);
    return M_PI*d*d/surfaceArea(a,c);
}

//- Projected area perpendicular to a direction making angle alpha
//  with the symmetry axis (exact for spheroids):
//  A_perp(alpha) = pi a sqrt(a^2 cos^2 + c^2 sin^2)
inline scalar projectedArea(const scalar a, const scalar c,
                            const scalar cosAlpha)
{
    const scalar c2 = cosAlpha*cosAlpha;
    return M_PI*a*std::sqrt(a*a*c2 + c*c*(1.0-c2));
}

//- Crosswise sphericity Phi_perp = (pi dEq^2/4)/A_perp(alpha)
inline scalar crosswiseSphericity(const scalar a, const scalar c,
                                  const scalar cosAlpha)
{
    const scalar d = dEq(a,c);
    return 0.25*M_PI*d*d/projectedArea(a,c,cosAlpha);
}

//- Lengthwise sphericity (Hoelzer-Sommerfeld):
//  Phi_par = (pi dEq^2 / 4) / (A_particle/2 - A_perp,mean_along)
//  For spheroids we use the standard approximation with the
//  projected area along the mean flow-parallel direction.
inline scalar lengthwiseSphericity(const scalar a, const scalar c,
                                   const scalar cosAlpha)
{
    const scalar d   = dEq(a,c);
    const scalar Apar = projectedArea(a,c,cosAlpha);      // area normal to flow
    const scalar half = 0.5*surfaceArea(a,c);
    const scalar denom = std::max(half - Apar, 1e-30);
    return 0.25*M_PI*d*d/denom;
}

//- Hoelzer & Sommerfeld (2008) orientation-dependent drag coefficient.
//  Re based on dEq and slip speed. cosAlpha = k.uhat.
inline scalar CdHolzerSommerfeld(const scalar Re,
                                 const scalar a, const scalar c,
                                 const scalar cosAlpha)
{
    const scalar Phi   = sphericity(a,c);
    const scalar PhiP  = crosswiseSphericity(a,c,cosAlpha);
    const scalar PhiL  = lengthwiseSphericity(a,c,cosAlpha);
    const scalar ReS   = std::max(Re, 1e-12);

    return 8.0/(ReS*std::sqrt(PhiL))
         + 16.0/(ReS*std::sqrt(Phi))
         + 3.0/(std::sqrt(ReS)*std::pow(Phi,0.75))
         + 0.42*std::pow(10.0, 0.4*std::pow(std::max(0.0,-std::log10(Phi)),0.2))/PhiP;
}

//- Di Felice voidage-correction exponent
inline scalar diFeliceBeta(const scalar Re)
{
    const scalar x = 1.5 - std::log10(std::max(Re,1e-12));
    return 3.7 - 0.65*std::exp(-0.5*x*x);
}

//- Bretherton constant B = (lam^2-1)/(lam^2+1), lam = c/a
inline scalar bretherton(const scalar a, const scalar c)
{
    const scalar l2 = (c/a)*(c/a);
    return (l2-1.0)/(l2+1.0);
}

//- Jeffery orbit period for validation case V3: T = 2 pi (lam + 1/lam)/gdot
inline scalar jefferyPeriod(const scalar a, const scalar c, const scalar gdot)
{
    const scalar lam = c/a;
    return 2.0*M_PI*(lam + 1.0/lam)/gdot;
}

// --- Lift / pitching torque coefficient shells -------------------------------
// Leading-order angular dependence ~ sin(a)cos(a); Re-dependent prefactors
// f(Re) must be calibrated per aspect ratio against DNS fits
// (Zastawny et al. 2012; Sanjeevi & Padding 2018). Values below are
// PLACEHOLDERS wired for the verification ladder, to be replaced by a small
// coefficient table read from the force-model dictionary.

inline scalar ClSpheroid(const scalar Re, const scalar sinA, const scalar cosA,
                         const scalar f1, const scalar f2)
{
    // C_L = (f1/Re^f2) sinA cosA   (dictionary-supplied f1,f2)
    return f1/std::pow(std::max(Re,1e-12), f2)*sinA*cosA;
}

inline scalar CtSpheroid(const scalar Re, const scalar sinA, const scalar cosA,
                         const scalar g1, const scalar g2)
{
    return g1/std::pow(std::max(Re,1e-12), g2)*sinA*cosA;
}


//- Rotational resistance functions for a prolate spheroid (Kim & Karrila).
//  Torque: T = 8 pi mu a^2 c [ XC dW_par + YC dW_perp + YH k x (E.k) ],
//  dW = 0.5 curl(u) - omega_p. Sphere limit: (XC,YC,YH) -> (1,1,0) so
//  T -> 8 pi mu R^3 (W_f - w_p). Identity YH/YC = Bretherton B(lam)
//  guarantees torque-free particles follow Jeffery's equation exactly.
inline void rotResist(const scalar a, const scalar c,
                      scalar& XC, scalar& YC, scalar& YH)
{
    if (c <= a*(1.0+1e-12)) { XC=1.0; YC=1.0; YH=0.0; return; }
    const scalar e = std::sqrt(1.0 - (a*a)/(c*c));
    const scalar L = std::log((1.0+e)/(1.0-e));
    XC = (4.0/3.0)*e*e*e*(1.0-e*e)/( 2.0*e - (1.0-e*e)*L);
    YC = (4.0/3.0)*e*e*e*(2.0-e*e)/(-2.0*e + (1.0+e*e)*L);
    YH = (4.0/3.0)*e*e*e*e*e      /(-2.0*e + (1.0+e*e)*L);
}
} // namespace sqGeom
} // namespace Foam

#endif
