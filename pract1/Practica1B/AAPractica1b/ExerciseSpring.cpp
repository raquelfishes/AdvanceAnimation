#include "Spring.h"
#include "Point.h"
#include "Scene.h"
#include "CGSolver.h"

///////////////////////////////////////////////////////////////////////////////////////////////////

void Spring::AddForces(bool implicit)
{
	//The compute of forces is the same in both methods (implicit and symplectic)
	Vec3 dLdX = (1.0f / length) * (a->GetPosition() - b->GetPosition());
	Vec3 force = (stiffness * (rest - length)) * dLdX;
	if (!HalfFixed()){
		//Both points not fixed
		a->AddForce(force);
		b->AddForce(-force);
	}
	else{
		//Some point fixed
		a->IsFixed() ? b->AddForce(-force) : a->AddForce(force);
	}
}

void Spring::AddToLinearSystem(CGSolver *solver, bool implicit)
{
	//Add fixed Matrix to A matrix 
	//A = M - h^2*K => K=dF/dp
	if (implicit){
		//If implicit method we have to compute rigid matrix
		//Compute dFdP from exercise 2 practice 1a
		//dFdp = (k*(L0 / L - 1)) * Matrix2::IDENTITY - ((k*L0 / L) * dLdx) * dLdx;
		Vec3 dLdX = (1.0f / length) * (a->GetPosition() - b->GetPosition());
		Matrix3 dFdP = (-1.0f) * pow(Scene::step, 2) * (stiffness * (rest / length - 1) * Matrix3::IDENTITY - (((stiffness * rest / length) * dLdX) ^ dLdX));
		
		solver->AddToDiagBlock(dFdP, a->GetSimIndex());
		solver->AddToDiagBlock(dFdP, b->GetSimIndex());
		solver->AddToSparseBlock(-dFdP, simIndex);
	}
}
