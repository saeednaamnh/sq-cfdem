// Unit tests of src/sqCfdem/common/spheroidGeometry.H against independent calculations.
#include <cstdio>
#include <random>
#include <algorithm>
#include "spheroidGeometry.H"
using namespace Foam; using namespace Foam::sqGeom;
static int fails=0, checks=0;
#define CHECK(cond,msg,...) do{checks++; if(!(cond)){fails++; printf("  FAIL: " msg "\n",__VA_ARGS__);}}while(0)
// independent references --------------------------------------------------
static double projNumeric(double a,double c,double ca,int n=1600){   // area of the projected ellipse, by direct sampling
  double sa=std::sqrt(std::max(0.0,1-ca*ca)); // view dir d = (sa,0,ca); axis along z
  // projection of spheroid onto plane normal to d: semi-axes a and sqrt(a^2 ca^2 + c^2 sa^2)
  // independent check: integrate the support-function form of the area of a convex body's shadow
  double sum=0; for(int i=0;i<n;i++){ double t=2*M_PI*(i+0.5)/n;
    // in-plane unit vectors: e1=(ca,0,-sa) , e2=(0,1,0); support h(u)=sqrt(a^2(ux^2+uy^2)+c^2 uz^2)
    double ux=std::cos(t)*ca, uy=std::sin(t), uz=-std::cos(t)*sa;
    double h=std::sqrt(a*a*(ux*ux+uy*uy)+c*c*uz*uz);
    // shadow area = 1/2 \oint (h^2 - h'^2) dt ; use numerical derivative of h
    double dt=1e-6, t2=t+dt; double ux2=std::cos(t2)*ca, uy2=std::sin(t2), uz2=-std::cos(t2)*sa;
    double h2=std::sqrt(a*a*(ux2*ux2+uy2*uy2)+c*c*uz2*uz2); double hp=(h2-h)/dt;
    sum+=0.5*(h*h-hp*hp)*(2*M_PI/n);} return sum; }
static double surfNumeric(double a,double c,int n=200000){ double s=0; // surface of revolution r(z)=a sqrt(1-z^2/c^2)
  for(int i=0;i<n;i++){ double th=M_PI*(i+0.5)/n; double r=a*std::sin(th), z=c*std::cos(th);
    double dr=a*std::cos(th), dz=-c*std::sin(th); s+=2*M_PI*r*std::sqrt(dr*dr+dz*dz)*(M_PI/n);} return s; }
static void kk(double lam,double&XC,double&YC,double&YH){ long double e=std::sqrt(1.0L-1.0L/(lam*(long double)lam)),
  L=std::log((1+e)/(1-e)); XC=(double)((4.0L/3)*e*e*e*(1-e*e)/(2*e-(1-e*e)*L));
  YC=(double)((4.0L/3)*e*e*e*(2-e*e)/(-2*e+(1+e*e)*L)); YH=(double)((4.0L/3)*e*e*e*e*e/(-2*e+(1+e*e)*L)); }
int main(){
  std::mt19937 g(7); std::uniform_real_distribution<double> U(-1,1), L(1.05,8.0), R(0.01,200);
  printf("[1] rotation matrix: orthonormal, det=+1, axis = R*z, for 2000 random unit quaternions\n");
  for(int i=0;i<2000;i++){ double w=U(g),x=U(g),y=U(g),z=U(g),n=std::sqrt(w*w+x*x+y*y+z*z); w/=n;x/=n;y/=n;z/=n;
    tensor M=quatToRot(w,x,y,z); double e=0;
    for(int r=0;r<3;r++)for(int c=0;c<3;c++){double d=0; for(int k=0;k<3;k++) d+=M(k,r)*M(k,c); e=std::max(e,std::fabs(d-(r==c)));}
    double det=M(0,0)*(M(1,1)*M(2,2)-M(1,2)*M(2,1))-M(0,1)*(M(1,0)*M(2,2)-M(1,2)*M(2,0))+M(0,2)*(M(1,0)*M(2,1)-M(1,1)*M(2,0));
    vector k=symmetryAxis(w,x,y,z); double ek=std::fabs(k.x()-M(0,2))+std::fabs(k.y()-M(1,2))+std::fabs(k.z()-M(2,2));
    CHECK(e<1e-12 && std::fabs(det-1)<1e-12 && ek<1e-14,"quaternion %d: orthonormality %.1e det %.15f axis %.1e",i,e,det,ek); }
  { double h=std::sqrt(0.5); vector k=symmetryAxis(h,-h,0,0); // LIGGGHTS "set quat 1 0 0 -90", as in the torque test
    CHECK(std::fabs(k.x())<1e-12 && std::fabs(k.y()-1)<1e-12 && std::fabs(k.z())<1e-12,"quat 1 0 0 -90 gives (%g,%g,%g), expected (0,1,0)",k.x(),k.y(),k.z());
    k=symmetryAxis(h,0,h,0); CHECK(std::fabs(k.x()-1)<1e-12&&std::fabs(k.z())<1e-12,"quat 0 1 0 90 gives (%g,%g,%g), expected (1,0,0)",k.x(),k.y(),k.z()); }
  printf("[2] projected area vs support-function integration, 300 random (lambda, angle)\n");
  double worst=0; for(int i=0;i<300;i++){ double lam=L(g),a=1,c=lam,ca=U(g); double r=projectedArea(a,c,ca)/projNumeric(a,c,ca);
    worst=std::max(worst,std::fabs(r-1)); } CHECK(worst<1e-6,"worst relative error %.2e",worst); printf("     worst relative error %.1e\n",worst);
  printf("[3] surface area vs numerical integration, lambda 1.001-20\n"); worst=0;
  for(double lam: {1.001,1.01,1.5,2.0,3.0,5.0,8.0,20.0}){ double r=surfaceArea(1,lam)/surfNumeric(1,lam); worst=std::max(worst,std::fabs(r-1)); }
  CHECK(worst<1e-8,"worst relative error %.2e",worst); printf("     worst relative error %.1e\n",worst);
  printf("[4] sphere limit: all sphericities 1, C_D = 24/Re + 3/sqrt(Re) + 0.42\n"); worst=0;
  for(double Re: {0.01,0.1,1.0,10.0,100.0,1000.0}) for(double ca: {1.0,0.5,0.0}){
    double cd=CdHolzerSommerfeld(Re,1,1,ca), ref=24/Re+3/std::sqrt(Re)+0.42; worst=std::max(worst,std::fabs(cd/ref-1)); }
  CHECK(worst<1e-5,"worst relative error %.2e",worst); printf("     worst relative error %.1e (exact sphere: Phi computed as %.17g)\n",worst,sphericity(1,1));
  printf("[5] resistance functions vs independent long-double Kim-Karrila, lambda 1.05-100; Bretherton identity\n"); worst=0; double wb=0;
  for(double lam=1.05; lam<100; lam*=1.07){ double XC,YC,YH,X2,Y2,H2; rotResist(1,lam,XC,YC,YH); kk(lam,X2,Y2,H2);
    worst=std::max({worst,std::fabs(XC/X2-1),std::fabs(YC/Y2-1),std::fabs(YH/H2-1)});
    wb=std::max(wb,std::fabs(YH/YC-bretherton(1,lam))); }
  CHECK(worst<1e-10,"worst relative error %.2e",worst); CHECK(wb<1e-12,"Bretherton identity error %.2e",wb);
  printf("     worst relative error %.1e, Bretherton identity error %.1e\n",worst,wb);
  printf("[6] near-sphere precision of the resistance functions (double vs long-double reference)\n");
  for(double d: {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6}){ double XC,YC,YH,X2,Y2,H2; rotResist(1,1+d,XC,YC,YH); kk(1+d,X2,Y2,H2);
    printf("     lambda = 1 + %.0e : X^C rel. err %.1e, Y^C %.1e  (values %.6f, %.6f)\n",d,std::fabs(XC/X2-1),std::fabs(YC/Y2-1),XC,YC); }
  printf("[7] Di Felice exponent: Re = 1 -> %.6f (reference 3.7 - 0.65 exp(-1.125) = %.6f)\n",diFeliceBeta(1.0),3.7-0.65*std::exp(-1.125));
  CHECK(std::fabs(diFeliceBeta(1.0)-(3.7-0.65*std::exp(-1.125)))<1e-12,"Di Felice %g",diFeliceBeta(1.0));
  printf("[8] oblate input (c < a): rotResist returns (X^C, Y^C, Y^H) = "); { double XC,YC,YH; rotResist(1,0.5,XC,YC,YH); printf("(%g, %g, %g)  <- sphere values, silently\n",XC,YC,YH); }
  printf("\n%d checks, %d failures\n",checks,fails); return fails!=0; }
