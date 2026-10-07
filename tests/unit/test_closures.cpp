#include <cstdio>
#include <random>
#include "spheroidGeometry.H"
using namespace Foam; using namespace Foam::sqGeom;
static int fails=0, checks=0;
#define CHECK(cond,...) do{checks++; if(!(cond)){fails++; printf("  FAIL: " __VA_ARGS__); printf("\n");}}while(0)
static std::vector<vector> oldMarkers(int nSat){ std::vector<vector> P; P.push_back(vector(0,0,0)); int nS=nSat-1; double g=M_PI*(3-std::sqrt(5.0));
  for(int i=0;i<nS;i++){ double r=std::cbrt((i+0.5)/nS), zz=1-2*(i+0.5)/nS, rho=std::sqrt(std::max(0.0,1-zz*zz)), th=g*i; P.push_back(vector(r*rho*std::cos(th),r*rho*std::sin(th),r*zz)); } return P; }
static double capFrac(double h){ return (2-3*h+h*h*h)/4; }   // exact volume fraction of the unit ball above z = h
int main(){
  const double rho=10, nu=1.5e-3, a=0.6785e-3, c=1.357e-3, U=0.05;
  printf("[1] drag without voidage reproduces the REM-verified V1 value; eps=1 with voidage equals no voidage\n");
  double Re, K=dragCoefficient(rho,nu,a,c,1.0,U,0.47,false,Re); printf("     F = %.4e N (V1 verified 1.2561e-05), Re = %.4f\n",K*U,Re);
  CHECK(std::fabs(K*U/1.2561e-5-1)<2e-4,"V1 force %.6e",K*U);
  double Re1, K1=dragCoefficient(rho,nu,a,c,0.3,U,1.0,true,Re1), Re0, K0=dragCoefficient(rho,nu,a,c,0.3,U,1.0,false,Re0);
  CHECK(std::fabs(K1/K0-1)<1e-14 && std::fabs(Re1-Re0)<1e-14,"eps=1: %.3e vs %.3e",K1,K0);
  printf("[2] with voidage: independent Di Felice formula, 500 random cases\n");
  std::mt19937 g(3); std::uniform_real_distribution<double> E(0.35,0.99),C(-1,1),L(1.0,6.0),V(1e-4,2.0);
  double worst=0; for(int i=0;i<500;i++){ double eps=E(g),ca=C(g),lam=L(g),u=V(g),aa=0.5e-3,cc=aa*lam; double Rv;
    double Kv=dragCoefficient(rho,nu,aa,cc,ca,u,eps,true,Rv);
    double d=2*std::cbrt(aa*aa*cc), Res=eps*d*u/nu, chi=3.7-0.65*std::exp(-0.5*std::pow(1.5-std::log10(Res),2));
    double Cd=CdHolzerSommerfeld(Res,aa,cc,ca), Kref=0.5*rho*Cd*0.25*M_PI*d*d*std::pow(eps,2-chi)*u;
    worst=std::max(worst,std::fabs(Kv/Kref-1)); }
  CHECK(worst<1e-13,"worst %.2e",worst); printf("     worst relative difference %.1e\n",worst);
  { double Rn,Rv; double Kn=dragCoefficient(rho,nu,a,c,1.0,U,1.0,false,Rn), Kv=dragCoefficient(rho,nu,a,c,1.0,U,0.47,true,Rv);
    printf("     at eps = 0.47 (your solid fraction 0.53): K_voidage / K_isolated = %.2f  (old code applied eps^-chi = %.1f)\n",Kv/Kn,std::pow(0.47,-diFeliceBeta(Rn))); }
  printf("[3] markers: centroid, second moments, radius-direction correlation, cap volumes\n");
  for(int n: {64,65}){ auto P=unitBallMarkers(n); double cx=0,cy=0,cz=0,xx=0,yy=0,zz=0,rmax=0,rz=0;
    for(auto&p:P){cx+=p.x();cy+=p.y();cz+=p.z();xx+=p.x()*p.x();yy+=p.y()*p.y();zz+=p.z()*p.z();rmax=std::max(rmax,std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z()));
      rz+=std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z())*p.z();}
    int N=P.size(); printf("     n=%d: count %d, centroid (%.1e, %.1e, %.1e), <xx>,<yy>,<zz> = %.3f, %.3f, %.3f, max radius %.3f, <r z> = %.1e\n",n,N,cx/N,cy/N,cz/N,xx/N,yy/N,zz/N,rmax,rz/N);
    CHECK(N==n && std::fabs(cx/N)<1e-15 && std::fabs(cy/N)<1e-15 && std::fabs(cz/N)<1e-15 && rmax<=1.0,"n=%d marker basics",n);
    CHECK(std::fabs(xx/N-0.2)<0.02 && std::fabs(yy/N-0.2)<0.02 && std::fabs(zz/N-0.2)<0.02,"n=%d moments",n); }
  { auto Pn=unitBallMarkers(64); auto Po=oldMarkers(64); double eN=0,eO=0;
    for(double h=-0.9; h<0.91; h+=0.05){ double fn=0,fo=0; for(auto&p:Pn) fn+=(p.z()>h); for(auto&p:Po) fo+=(p.z()>h);
      eN=std::max(eN,std::fabs(fn/Pn.size()-capFrac(h))); eO=std::max(eO,std::fabs(fo/Po.size()-capFrac(h))); }
    printf("     volume above a plane at height h along the axis, worst error over h: new %.3f, old %.3f\n",eN,eO);
    CHECK(eN<eO,"new markers not better than old: %.3f vs %.3f",eN,eO); }
  printf("[4] shape check\n");
  struct T{double a,b,c,n1,n2; bool ok;} tc[]={{1,1,3,2,2,true},{1,1,1,2,2,true},{1,1,0.5,2,2,false},{1,1.2,3,2,2,false},{1,1,3,4,2,false},{1,1,3,2,8,false}};
  for(auto&t:tc){ std::string why=unsupportedShape(t.a,t.b,t.c,t.n1,t.n2); bool ok=why.empty();
    printf("     a=%g b=%g c=%g n=(%g,%g): %s%s\n",t.a,t.b,t.c,t.n1,t.n2,ok?"accepted":"rejected: ",why.c_str()); CHECK(ok==t.ok,"shape case"); }
  printf("\n%d checks, %d failures\n",checks,fails); return fails!=0; }
