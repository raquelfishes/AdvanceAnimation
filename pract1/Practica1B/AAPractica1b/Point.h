#pragma once

#include <Inventor/SbLinear.h>
#include "Vec3.h"

class CGSolver;

class Point
{
public:
   static float mass;

protected:
   Vec3 position;
   Vec3 velocity;

   Vec3 force;

   bool fixed;
   int simIndex;

public:
   Point(void);
   ~Point(void);

   void InitPositionCorner(int i, int j, int k, int x, int y, int z);
   void InitPositionCenter(int i, int j, int k, int x, int y, int z);

   void ClearForce(void);
   void AddGravity(void);
   void AddForce(const Vec3& f);
   void GetState(float *x, float *v) const;
   void SetState(float *x, float *v);
   const Vec3& GetPosition(void) const
   {
      return position;
   }
   SbVec3f GetPositionRender(void) const
   {
      SbVec3f pos(position[0], position[1], position[2]);
      return pos;
   }
   const Vec3& GetVelocity(void) const
   {
      return velocity;
   }
   const Vec3& GetForce(void) const
   {
	   return force;
   }
   int GetSimIndex(void) const
   {
      return simIndex;
   }
   void SetSimIndex(int ind)
   {
      simIndex = ind;
   }
   bool IsFixed(void) const
   {
      return fixed;
   }
   void SetFixed(bool f)
   {
      fixed = f;
   }

   void AddToLinearSystem(CGSolver *solver, bool implicit);

   void Update(void);

};
