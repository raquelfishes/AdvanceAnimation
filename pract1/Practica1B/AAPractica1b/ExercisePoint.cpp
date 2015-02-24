#include "Point.h"
#include "Scene.h"
#include "CGSolver.h"

//Method to add a value on the diagonal of the matrix
void addDiagonal(float value, Matrix3& matrix){
	matrix[0][0] = value;
	matrix[1][1] = value;
	matrix[2][2] = value;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void Point::AddToLinearSystem(CGSolver *solver, bool implicit)
{
	//The compute of A matrix and b vector is the same in both methods (implicit and symplectic)
	//A*v = b
	//A = M
	//b = M*v+h*F
	//Adding Mass matrix block
	Matrix3 mass = Matrix3::ZERO;
	addDiagonal(this->mass, mass);
	solver->AddToDiagBlock(mass, simIndex);
	//Adding vector b
	Vec3 b = Vec3::ZERO;
	b = mass*velocity;
	b += Scene::step*force;
	solver->AddToVectorBlock(b, simIndex);
}

