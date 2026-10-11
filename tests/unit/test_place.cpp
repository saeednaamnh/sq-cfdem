#include <cstdio>
#include <random>
#include "spheroidGeometry.H"
using namespace Foam::sqGeom;
int main(){ int fails=0, n=0; const double lo[3]={0,0,0}, hi[3]={15,10,10}; const bool per[3]={true,true,false};
  auto T=[&](double x,double y,double z,int want,double ex,double ey){ double p[3]={x,y,z}; int r=placeMarker(p,lo,hi,per); n++;
    bool ok = r==want && (r==2 || (std::fabs(p[0]-ex)<1e-12 && std::fabs(p[1]-ey)<1e-12 && p[2]==z));
    if(!ok){fails++; printf("  FAIL (%g,%g,%g): got %d -> (%g,%g,%g), want %d (%g,%g)\n",x,y,z,r,p[0],p[1],p[2],want,ex,ey);} };
  T(7,5,5, 0, 7,5);            // inside
  T(-0.3,5,5, 1, 14.7,5);      // beyond x = 0 -> wrapped to the far side
  T(15.2,5,5, 1, 0.2,5);       // beyond x = Lx
  T(15,5,5, 1, 0,5);           // exactly on the upper periodic face -> lower face
  T(3,-0.5,5, 1, 3,9.5);       // beyond y = 0
  T(-0.3,10.4,5, 1, 14.7,0.4); // corner: both wrapped
  T(7,5,-0.1, 2, 0,0);         // below the bottom wall -> outside the fluid
  T(7,5,10.2, 2, 0,0);         // above the top wall
  T(-31,25,5, 1, 14,5);        // several periods away
  std::mt19937 g(5); std::uniform_real_distribution<double> U(-20,40), Z(0,10); int bad=0;
  for(int i=0;i<100000;i++){ double p[3]={U(g),U(g),Z(g)}; placeMarker(p,lo,hi,per); if(p[0]<0||p[0]>=15||p[1]<0||p[1]>=10) bad++; }
  n++; if(bad){fails++; printf("  FAIL: %d random points not wrapped into the box\n",bad);}
  printf("%d checks, %d failures\n",n,fails); return fails!=0; }
