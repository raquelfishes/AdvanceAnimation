#include <iostream>
using namespace std;
#include "Exercise.h"
#include "LinSys.h"

// gravitational acceleration (9.81)
static const float g = 9.81f;

// spring stiffness
static const float k = 100.0f;

// damping
static const float d = 0.01f;

void AdvanceTimeStep(vector<RigidBody>& bodies, float step, bool constraint, bool collisions, float floor)
{
}

// Linear system example
void LinearSystemExample(void)
{
   //Set up matrix A and vector b
   int size=3;
   MatrixMN A(size, size);
   A(0,0)=1.0f; A(0,1)=2.0f; A(0,2)=3.0f;
   A(1,0)=4.0f; A(1,1)=5.0f; A(1,2)=6.0f;
   A(2,0)=7.0f; A(2,1)=8.0f; A(2,2)=9.0f;
   Vector b(size);
   b(0)=10.0f; b(1)=11.0f; b(2)=12.0f;
   Vector x(size);

   //Solve linear system A*x=b
   GaussElimination(A, b, x);
}

// LCP example
void LCPExample(void)
{
   //Set up matrix A and vector b
   int size=3;
   MatrixMN A(size, size);
   A(0,0)=1.0f; A(0,1)=2.0f; A(0,2)=3.0f;
   A(1,0)=4.0f; A(1,1)=5.0f; A(1,2)=6.0f;
   A(2,0)=7.0f; A(2,1)=8.0f; A(2,2)=9.0f;
   Vector b(size);
   b(0)=10.0f; b(1)=11.0f; b(2)=12.0f;
   Vector x(size);

   //Set up vector of flags. flags[i] = true indicates inequality constrain; flags[i] = false indicates equality constraint
   bool *flags = new bool[size];
   flags[0] = false; flags[1] = true; flags[2] = false;

   //Solve linear system A*x >= b, x >= 0, x^T (A*x - b) = 0
   ProjectedGaussSeidel(A, b, x, flags);
}




