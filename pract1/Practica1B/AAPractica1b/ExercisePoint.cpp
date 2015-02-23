#include "Point.h"
#include "Scene.h"
#include "CGSolver.h"

void addImplicitSystem(Point* point, CGSolver *solver){

}

void addDiagonal(float value, Matrix3 matrix){
	matrix[0][0] = value;
	matrix[1][1] = value;
	matrix[2][2] = value;
}

void addSymplecticSystem(Point* point, CGSolver *solver){
	//Adding Mass matrix block
	Matrix3 mass = Matrix3::ZERO;
	addDiagonal(point->mass,mass);
	solver->AddToDiagBlock(mass, point->GetSimIndex());
	//Adding vector b
	Vec3 b = Vec3::ZERO;
	b = mass*point->GetVelocity();
	b += Scene::step*point->GetForce();
	solver->AddToVectorBlock(b, point->GetSimIndex());

}

///////////////////////////////////////////////////////////////////////////////////////////////////

void Point::AddToLinearSystem(CGSolver *solver, bool implicit)
{
	if (implicit){
		//Calculate implicit forces
		addImplicitSystem(this,solver);
	}
	else{
		//Calculate symplectic forces
		addSymplecticSystem(this,solver);
	}
}

