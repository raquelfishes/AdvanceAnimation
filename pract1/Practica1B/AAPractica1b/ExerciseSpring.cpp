#include "Spring.h"
#include "Point.h"
#include "Scene.h"
#include "CGSolver.h"

///////////////////////////////////////////////////////////////////////////////////////////////////

Vec3 forceSpring(Spring* spring){
	//Compute force spring from exercise 2 practice 1b
	//L = (p1 - p2).length();
	//dLdx = (1.0f / L) * (p1 - p2);
	//F = (k*(L0 - L)) * dLdx;
	Vec3 dLdX = (1.0f / spring->GetLength()) * (spring->GetPointA()->GetPosition() - spring->GetPointB()->GetPosition());
	Vec3 force = (spring->stiffness * (spring->GetRest() - spring->GetLength())) * dLdX;
	return force;
}


void addSymplecticForces(Spring* spring){
	if (!spring->HalfFixed()){
		//Both points not fixed
		Vec3 force = forceSpring(spring);
		spring->GetPointA()->AddForce(force);
		spring->GetPointB()->AddForce(-force);
	}
	else{
		//Some point fixed
		Vec3 force = forceSpring(spring);
		if (spring->GetPointA()->IsFixed()){
			spring->GetPointB()->AddForce(-force);
		}
		else{
			spring->GetPointA()->AddForce(force);
		}
	}
}

void addImplicitForces(Spring* spring){
	Vec3 force = forceSpring(spring);
	if (!spring->HalfFixed()){
		//Both points not fixed
		spring->GetPointA()->AddForce(force);
		spring->GetPointB()->AddForce(-force);
	}
	else{
		//Some point fixed
		if (spring->GetPointA()->IsFixed()){
			spring->GetPointB()->AddForce(-force);
		}
		else{
			spring->GetPointA()->AddForce(force);
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////

void Spring::AddForces(bool implicit)
{
	if (implicit){
		//Calculate implicit forces
		addImplicitForces(this);
	}
	else{
		//Calculate symplectic forces
		addSymplecticForces(this);
	}
}

void Spring::AddToLinearSystem(CGSolver *solver, bool implicit)
{
}
