#include "GrabForce.h"
#include "Point.h"
#include "Spring.h"
#include "CGSolver.h"
#include "Scene.h"


GrabForce::GrabForce(void)
{
}

GrabForce::GrabForce(Point *p)
{
   point=p;
   grabPos=point->GetPosition();
}

GrabForce::~GrabForce(void)
{
}

void GrabForce::UpdateGrab(const Vec3& pos)
{
   grabPos=pos;
}

void GrabForce::ApplyForce(void)
{
   //Force proportional to the distance between the current mouse position and the position of the point
   point->AddForce(Spring::stiffness*(grabPos-point->GetPosition()));
}

void GrabForce::AddToLinearSystem(CGSolver *solver)
{
   //Compute value
   float val=Spring::stiffness*Scene::step*Scene::step;

   //Add values in the diagonal band of A
   solver->AddToDiagBlock(val*Matrix3::IDENTITY, point->GetSimIndex());
 
}

