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

//axis with distance to node
Vector3 axis(-0.5, 0, 0);

//Similar a Spring2D from practice 1
class Spring{
	
private:
	const Vector3 p1;
	const Vector3 p2;
	Vector3 f;
	const float k;

public:
	Spring(const Vector3 p1, const Vector3 p2, const float k) :p1(p1), p2(p2), k(k){ f = -k*(p1 - p2); }
	~Spring(void){}
	const Vector3 getF1(void) const { return f; }
	Vector3 getF2(void) const { return -f; }
};

Matrix33 fixRotation(const Matrix33 rot){
	
	Vector3 alfa = rot.getColumn(0);
	alfa.normalize();

	Vector3 beta = rot.getColumn(1);
	beta -= beta.dotProd(beta, alfa)*alfa;

	Vector3 gamma = Vector3::crossProd(alfa, beta);

	return Matrix33(alfa, beta, gamma);
}


void advanceConstraints(vector<RigidBody>& bodies, float step, bool collisions, float floor){

}

void advanceSprings(vector<RigidBody>& bodies, float step, bool collisions, float floor){
	int nbodies = bodies.size();

	//Initialize forces with gravity and damping
	vector<Vector3> F(nbodies), T(nbodies);
	for (int i = 0; i < nbodies; i++){
		F[i] = Vector3(0.0f, -bodies[i].Mass()*g, 0.0f);
		F[i] -= d*bodies[i].Velocity();
		T[i] = -d*bodies[i].Omega();
	}

	//Add springs forces
	//Create first spring with position 0 as p1
	Spring sp0 (Vector3::ZERO, bodies[0].posLocalToGlobal(axis), k);
	F[0] += sp0.getF2();
	T[0] += Vector3::crossProd(bodies[0].vecLocalToGlobal(axis),sp0.getF2());

	//Add spring force for the rest of bodies
	for (int i = 1; i < nbodies; i++){
		Spring sp(bodies[i - 1].posLocalToGlobal(-axis), bodies[i].posLocalToGlobal(axis), k);
		F[i - 1] += sp.getF1();
		F[i] += sp.getF2();
		T[i - 1] += Vector3::crossProd(bodies[i - 1].vecLocalToGlobal(-axis), sp.getF1());
		T[i] += Vector3::crossProd(bodies[i].vecLocalToGlobal(axis), sp.getF2());
	}

	//Collisions
	if (collisions){
		for (int i = 0; i < nbodies; i++){
			//Check first node of body
			Vector3 positionA = bodies[i].posLocalToGlobal(axis);
			if (positionA[1] < floor){
				//Add force floor
				Vector3 fFloor = (k*(floor - positionA[1]))*Vector3::UNIT_Y;
				F[i] += fFloor;
				T[i] += Vector3::crossProd(bodies[i].vecLocalToGlobal(axis), fFloor);
			}
			//Check second node of body
			Vector3 positionB = bodies[i].posLocalToGlobal(-axis);
			if (positionB[1] < floor){
				//Add force floor
				Vector3 fFloor = (k*(floor - positionB[1]))*Vector3::UNIT_Y;
				F[i] += fFloor;
				T[i] += Vector3::crossProd(bodies[i].vecLocalToGlobal(-axis), fFloor);
			}
		}
	}

	//Integrate
	//velocities
	for (int i = 0; i < nbodies; i++){
		// v <- v + dt/m * F
		bodies[i].SetVelocity(bodies[i].Velocity() + (step / bodies[i].Mass())*F[i]);
		// w <- w + dt * M^-1 * (T - w x Mw )
		bodies[i].SetOmega(bodies[i].Omega() + step*bodies[i].InertiaInv()*(T[i] + ((-1.0f)*Vector3::crossProd(bodies[i].Omega(), bodies[i].Inertia()*bodies[i].Omega()))));
	}

	//positions
	for (int i = 0; i < nbodies; i++){
		// x <- x + dt * v
		bodies[i].SetPosition(bodies[i].Position() + step*bodies[i].Velocity());
		// R <- R + dt*w*R
		bodies[i].SetRotation(fixRotation(bodies[i].Rotation() + step*Matrix33::ToCrossProd(bodies[i].Omega())*bodies[i].Rotation()));
		bodies[i].UpdateInertia();
	}
}


void AdvanceTimeStep(vector<RigidBody>& bodies, float step, bool constraint, bool collisions, float floor)
{
	(constraint) ? advanceConstraints(bodies, step, collisions, floor) : advanceSprings(bodies, step, collisions, floor);
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




