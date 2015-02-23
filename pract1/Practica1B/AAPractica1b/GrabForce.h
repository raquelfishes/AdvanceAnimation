#pragma once

#include "Vec3.h"

//Forward declarations
class Point;
class CGSolver;

class GrabForce
{
public:
   Point *point;

   Vec3 grabPos;

public:
   GrabForce(void);
   GrabForce(Point *p);
   ~GrabForce(void);

   void UpdateGrab(const Vec3& pos);

   void ApplyForce(void);
   void AddToLinearSystem(CGSolver *solver);
};
