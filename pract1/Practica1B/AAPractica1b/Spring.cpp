#include "Spring.h"
#include "Point.h"
#include "Scene.h"
#include "CGSolver.h"


float Spring::stiffness=1.0f;

Spring::Spring(void)
{
}

Spring::~Spring(void)
{
}

void Spring::Initialize(Point *a_, Point *b_)
{
   //set end points
   a=a_;
   b=b_;

   //compute rest length
   rest=ComputeLength();

   //compute mid pos
   ComputePosition();
}

float Spring::ComputeLength(void)
{
   length=(a->GetPosition()-b->GetPosition()).length();
   return length;
}

void Spring::ComputePosition(void)
{
   position=0.5f*a->GetPositionRender()+0.5f*b->GetPositionRender();
}

void Spring::Update(void)
{
   ComputeLength();
   ComputePosition();
}

bool Spring::FullFixed(void)
{
   return a->IsFixed() && b->IsFixed();
}

bool Spring::HalfFixed(void)
{
   return a->IsFixed() || b->IsFixed();
}

void Spring::InitLinearSystem(CGSolver *solver)
{
   if(!a->IsFixed() && !b->IsFixed())
   {
      solver->InitSparseBlockIndices(simIndex, a->GetSimIndex(), b->GetSimIndex());
   }
}