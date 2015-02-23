#pragma once

#include <Inventor/SbLinear.h>

//Forward declarations
class Point;
class CGSolver;

class Spring
{
public:
   static float stiffness;

protected:
   Point *a, *b;

   float rest;
   float length;

   SbVec3f position;

   int simIndex;

public:
   Spring(void);
   ~Spring(void);

   void Initialize(Point *a_, Point *b_);
   float ComputeLength(void);
   void ComputePosition(void);
   const SbVec3f& GetPosition(void) const
   {
      return position;
   }
   float GetLength(void) const
   {
      return length;
   }
   float GetRest(void) const
   {
	   return rest;
   }
   Point *GetPointA(void) const
   {
      return a;
   }
   Point *GetPointB(void) const
   {
      return b;
   }
   int GetSimIndex(void) const
   {
      return simIndex;
   }
   void SetSimIndex(int ind)
   {
      simIndex = ind;
   }

   void AddForces(bool implicit);

   void AddToLinearSystem(CGSolver *solver, bool implicit);

   void InitLinearSystem(CGSolver *solver);

   void Update(void);

   bool HalfFixed(void);
   bool FullFixed(void);
};
