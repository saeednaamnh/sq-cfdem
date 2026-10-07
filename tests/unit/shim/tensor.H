#pragma once
#include "vector.H"
namespace Foam { struct tensor { scalar t[9];
  tensor(scalar a,scalar b,scalar c,scalar d,scalar e,scalar f,scalar g,scalar h,scalar i){
   t[0]=a;t[1]=b;t[2]=c;t[3]=d;t[4]=e;t[5]=f;t[6]=g;t[7]=h;t[8]=i;}
  scalar operator()(int r,int c) const {return t[3*r+c];} }; }
