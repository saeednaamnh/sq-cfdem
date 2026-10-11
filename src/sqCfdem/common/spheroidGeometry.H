/*---------------------------------------------------------------------------*\
    sqCfdem - superquadric extension for CFDEMcoupling-PUBLIC
    License: GPLv3 (inherited from CFDEMcoupling)

    spheroidGeometry.H
    Header-only helpers for prolate spheroids (a = b <= c), n1 = n2 = 2.
    All functions are pure math: no OpenFOAM state, unit-agnostic.
    Models must call unsupportedShape() and abort on a non-empty reason:
    the closures below are valid for prolate spheroids and spheres only.
\*---------------------------------------------------------------------------*/

#ifndef spheroidGeometry_H
#define spheroidGeometry_H

#include "vector.H"
#include "tensor.H"
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

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

//- Mean longitudinal projected area (Hoelzer & Sommerfeld 2008): the projection
//  of the particle onto planes PARALLEL to the flow, averaged over all viewing
//  directions normal to the flow. With the symmetry axis at angle alpha to the
//  flow, a viewing direction at azimuth theta makes cos(beta) = sin(alpha) cos(theta)
//  with the axis, and the projection is projectedArea(a, c, cos(beta)).
inline scalar meanLongitudinalArea(const scalar a, const scalar c,
                                   const scalar cosAlpha)
{
    const scalar sinA = std::sqrt(std::max(0.0, 1.0 - cosAlpha*cosAlpha));
    const int n = 64;                      // midpoint rule on theta in [0, pi/2]
    scalar sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const scalar th = (i + 0.5)*0.5*M_PI/n;
        sum += projectedArea(a, c, sinA*std::cos(th));
    }
    return sum/n;
}

//- Lengthwise sphericity (Hoelzer & Sommerfeld 2008):
//  Phi_par = (pi dEq^2/4) / (A_particle/2 - mean longitudinal projected area)
inline scalar lengthwiseSphericity(const scalar a, const scalar c,
                                   const scalar cosAlpha)
{
    const scalar d     = dEq(a,c);
    const scalar Along = meanLongitudinalArea(a,c,cosAlpha);
    const scalar half  = 0.5*surfaceArea(a,c);
    const scalar denom = std::max(half - Along, 1e-30);
    return 0.25*M_PI*d*d/denom;
}

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

//- Supported shapes: superquadrics with blockiness n1 = n2 = 2 and semi-axes
//  a = b <= c (prolate spheroids, spheres). Returns an empty string if the shape
//  is supported, otherwise the reason; callers abort with it, because the drag,
//  torque and void-fraction closures are not valid for other shapes.
inline std::string unsupportedShape(const scalar a, const scalar b, const scalar c,
                                    const scalar n1, const scalar n2)
{
    const scalar tol = 1e-6;
    if (std::fabs(n1 - 2.0) > tol || std::fabs(n2 - 2.0) > tol)
        return "blockiness (n1, n2) must be (2, 2): the drag and torque closures "
               "are for spheroids, not general superquadrics";
    if (std::fabs(a - b) > tol*std::max(a, b))
        return "semi-axes a and b must be equal: the closures are for "
               "axisymmetric spheroids";
    if (c < a*(1.0 - tol))
        return "oblate spheroids (c < a) are not supported: the closures are "
               "for prolate spheroids and spheres";
    return std::string();
}

//- Drag coefficient K such that F = K (u_f - u_p), with the Hoelzer &
//  Sommerfeld (2008) C_D referred to the volume-equivalent sphere.
//  Without voidage correction: Re = d |u_r| / nu, K = 1/2 rho C_D A |u_r|.
//  With the Di Felice (1994) correction the single-particle drag is evaluated
//  with the SUPERFICIAL slip eps |u_r| (Reynolds number included) and multiplied
//  by eps^-chi, so K = 1/2 rho C_D(Re) A eps^(2-chi) |u_r| with Re = eps d |u_r|/nu,
//  exactly as in CFDEMcoupling's DiFeliceDrag. Re returns the value used.
inline scalar dragCoefficient
(
    const scalar rho, const scalar nu, const scalar a, const scalar c,
    const scalar cosAlpha, const scalar magUr, const scalar eps,
    const bool voidage, scalar& Re
)
{
    const scalar d    = dEq(a,c);
    const scalar Aref = 0.25*M_PI*d*d;
    if (!voidage)
    {
        Re = d*magUr/nu;
        return 0.5*rho*CdHolzerSommerfeld(Re, a, c, cosAlpha)*Aref*magUr;
    }
    Re = eps*d*magUr/nu;
    return 0.5*rho*CdHolzerSommerfeld(Re, a, c, cosAlpha)*Aref*magUr
          *std::pow(eps, 2.0 - diFeliceBeta(Re));
}

//- Satellite markers filling the unit ball, for the divided void fraction.
//  Antipodal pairs, so the centroid is exactly zero and the deposited volume
//  stays centred on the particle for any orientation; directions on a
//  Fibonacci hemisphere, radii from equal-volume shells assigned in an order
//  decorrelated from the direction. An odd count adds the centre point.
inline std::vector<vector> unitBallMarkers(const int n)
{
    std::vector<vector> P;
    if (n % 2) P.push_back(vector(0, 0, 0));
    const int m = n/2;                                 // antipodal pairs
    if (m == 0) return P;
    int step = std::max(1, int(std::floor(0.6180339887*m)));
    auto gcd = [](int x, int y) { while (y) { const int t = x % y; x = y; y = t; } return x; };
    while (gcd(step, m) != 1) ++step;                  // permutation j -> j*step mod m
    const scalar golden = M_PI*(3.0 - std::sqrt(5.0));
    for (int j = 0; j < m; j++)
    {
        const scalar zz  = 1.0 - (j + 0.5)/m;          // upper hemisphere; antipode covers the lower
        const scalar rho = std::sqrt(std::max(0.0, 1.0 - zz*zz));
        const scalar th  = golden*j;
        const scalar r   = std::cbrt((((j*step) % m) + 0.5)/m);
        const vector p(r*rho*std::cos(th), r*rho*std::sin(th), r*zz);
        P.push_back(p);
        P.push_back(vector(-p.x(), -p.y(), -p.z()));
    }
    return P;
}

//- Place a marker relative to the fluid domain [lo, hi]. In a periodic
//  direction a point outside the box is wrapped back into it (x is modified);
//  in a non-periodic direction (walls) a point outside the box lies outside the
//  fluid, and its share of particle volume belongs to no cell.
//  Returns 0: inside, unchanged; 1: inside after wrapping; 2: outside the domain.
enum markerPlace { markerInside = 0, markerWrapped = 1, markerOutside = 2 };
inline int placeMarker(scalar x[3], const scalar lo[3], const scalar hi[3],
                       const bool periodic[3])
{
    bool wrapped = false;
    for (int d = 0; d < 3; d++)
    {
        const scalar L = hi[d] - lo[d];
        if (periodic[d])
        {
            if (x[d] < lo[d] || x[d] >= hi[d])
            {
                const scalar f = (x[d] - lo[d])/L;
                x[d] = lo[d] + L*(f - std::floor(f));
                if (x[d] >= hi[d]) x[d] = lo[d];          // guard against round-off
                wrapped = true;
            }
        }
        else if (x[d] < lo[d] || x[d] > hi[d])
        {
            return markerOutside;
        }
    }
    return wrapped ? markerWrapped : markerInside;
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


//- Rotational resistance functions for a prolate spheroid (Kim & Karrila),
//  normalised by the semi-major axis c.
//  Torque: T = 8 pi mu c^3 [ XC dW_par + YC dW_perp + YH k x (E.k) ],
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
