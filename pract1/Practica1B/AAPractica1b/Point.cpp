#include "Point.h"
#include "Scene.h"
#include "CGSolver.h"

#include <iostream>
using namespace std;


float Point::mass=1.0f;

Point::Point(void)
{
   fixed=false;

   velocity = Vec3::ZERO;
}

Point::~Point(void)
{
}

void Point::InitPositionCorner(int i, int j, int k, int x, int y, int z)
{
   position = Vec3(Scene::xSize*(float)i/(float)(x-1),
      Scene::ySize*(float)j/(float)(y-1), Scene::zSize*(float)k/(float)(z-1));
}

void Point::InitPositionCenter(int i, int j, int k, int x, int y, int z)
{
   float px=Scene::xSize*(float)(2*i+1)/(float)(2*x-2);
   float py=Scene::ySize*(float)(2*j+1)/(float)(2*y-2);
   float pz=Scene::zSize*(float)(2*k+1)/(float)(2*z-2);
   position = Vec3(px, py, pz);
}

void Point::ClearForce(void)
{
   force = Vec3::ZERO;
}

void Point::AddGravity(void)
{
   force+=Vec3(0.0f, -mass*9.8f, 0.0f);
}

void Point::AddForce(const Vec3& f)
{
   force+=f;
}

void Point::GetState(float *x, float *v) const
{
   x[0] = position[0]; x[1] = position[1]; x[2] = position[2];
   v[0] = velocity[0]; v[1] = velocity[1]; v[2] = velocity[2];
}

void Point::SetState(float *x, float *v)
{
   position = Vec3(x[0], x[1], x[2]);
   velocity = Vec3(v[0], v[1], v[2]);
}

void Point::Update(void)
{
   if(fixed)
   {
      return;
   }
   cout << "valores antes:" << velocity[0] << " " << velocity[1] << " " << velocity[2] << "//" << position[0] << " " << position[1] << " " << position[2] << endl;
   velocity+=(1.0f/mass*Scene::step)*force;
   position+=Scene::step*velocity;
   cout << "valores despues:" << velocity[0] << " " << velocity[1] << " " << velocity[2] << "//" << position[0] << " " << position[1] << " " << position[2] << endl;

}