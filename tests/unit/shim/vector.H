#pragma once
#include <cmath>
namespace Foam { typedef double scalar;
struct vector { scalar v[3]; vector(){} vector(scalar a,scalar b,scalar c){v[0]=a;v[1]=b;v[2]=c;}
  scalar x() const{return v[0];} scalar y() const{return v[1];} scalar z() const{return v[2];} }; }
