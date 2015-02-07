#pragma once

#include <Inventor/SbLinear.h>

class Sphere
{
public:
   SbVec3f position;

   bool fixed;

public:
   Sphere(void){fixed=false;}
   ~Sphere(void){}

};


class Tube
{
public:
   Sphere *x, *y;

   float length;
   SbVec3f position;

public:
   Tube(void){}
   ~Tube(void){}

   void Initialize(Sphere *a, Sphere *b)
   {
      x=a;
      y=b;
      Update();
   }

   void Update(void)
   {
      length=(x->position-y->position).length();
      position=0.5f*x->position+0.5f*y->position;
   }

};

