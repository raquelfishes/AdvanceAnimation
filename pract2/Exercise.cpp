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

//Class to compute the joint constraint
class JointConstraint
{
private:
	const RigidBody& b1, b2;
	int bodyId1, bodyId2;
	const Vector3& r1, r2;
	Vector3 C, v;
	Matrix33 Jw1, Jw2;

public:
	JointConstraint(void) : b1(RigidBody::RIGIDBODYVOID), b2(RigidBody::RIGIDBODYVOID), r1(Vector3::ZERO), r2(Vector3::ZERO) {}
	JointConstraint(const RigidBody& b1, const RigidBody& b2, int bodyId1, int bodyId2, const Vector3& r1, const Vector3& r2)
		: b1(b1), b2(b2), bodyId1(bodyId1), bodyId2(bodyId2), r1(r1), r2(r2)
	{
		// Restriction C => b1 = b2     => They must be the same point
		C = b1.posLocalToGlobal(r1) - b2.posLocalToGlobal(r2);
		v = b1.pointVelocity(r1) - b2.pointVelocity(r2);
		Jw1 = Matrix33::ToCrossProd(b1.vecLocalToGlobal(-r1));
		Jw2 = Matrix33::ToCrossProd(b2.vecLocalToGlobal(r2));
	}
	int GetBodyId1(void) const { return bodyId1; }
	int GetBodyId2(void) const { return bodyId2; }
	const Vector3& GetC(void) const { return C; }
	const Vector3& GetV(void) const { return v; }
	const Matrix33& GetJv1(void) const { return Matrix33::IDENTITY; }
	Matrix33 GetJv2(void) const { return -Matrix33::IDENTITY; }
	const Matrix33& GetJw1(void) const { return Jw1; }
	const Matrix33& GetJw2(void) const { return Jw2; }
};

//Class to compute the joint constraint at the first body. This body is fixed!
class JointFixedConstraint
{
private:
	const RigidBody& b;
	int bodyId;
	const Vector3& r, p;
	Vector3 C, v;
	Matrix33 Jw;

public:
	JointFixedConstraint(void) : b(RigidBody::RIGIDBODYVOID), r(Vector3::ZERO), p(Vector3::ZERO) {}
	JointFixedConstraint(const RigidBody& b, int bodyId, const Vector3& r, const Vector3& p) : b(b), bodyId(bodyId), r(r), p(p)
	{
		// Restriction C => p(0,0,0) = b     => They must be the same point
		C = p - b.posLocalToGlobal(r);
		v = -b.pointVelocity(r);
		Jw = Matrix33::ToCrossProd(b.vecLocalToGlobal(r));
	}
	int GetBodyId(void) const { return bodyId; }
	const Vector3& GetC(void) const { return C; }
	const Vector3& GetV(void) const { return v; }
	Matrix33 GetJv(void) const { return -Matrix33::IDENTITY; }
	const Matrix33& GetJw(void) const { return Jw; }
};

//Class to compute contact to floor constraint
class ContactConstraint
{
private:
	const RigidBody& b;
	const Vector3& r;
	int bodyId;
	float p;
	float C, v;
	Vector3 Jw;

public:
	ContactConstraint(void) : b(RigidBody::RIGIDBODYVOID), r(Vector3::ZERO) {}
	ContactConstraint(const RigidBody& b, int bodyId, const Vector3& r, float p)
		: b(b), bodyId(bodyId), r(r), p(p)
	{
		// Restriction C => p(0,0,0) = b     => They must be the same point
		C = b.posLocalToGlobal(r)[1] - p;
		v = b.pointVelocity(r)[1];
		Jw = Matrix33::ToCrossProd(b.vecLocalToGlobal(r)) * Vector3::UNIT_Y;
	}
	int GetBodyId(void) const { return bodyId; }
	const float& GetC(void) const { return C; }
	const float& GetV(void) const { return v; }
	const Vector3& GetJv(void) const { return Vector3::UNIT_Y; }
	const Vector3& GetJw(void) const { return Jw; }
};

// Orthonormalization or rotation
Matrix33 fixRotation(const Matrix33 rot){
	// R = (alfa beta gamma)
	// alfa <- alfa / ||alfa||
	Vector3 alfa = rot.getColumn(0);
	alfa.normalize();
	// beta <- beta - (beta^T*alfa)*alfa
	Vector3 beta = rot.getColumn(1);
	beta -= beta.dotProd(beta, alfa)*alfa;
	// gamma <- alfa x beta
	Vector3 gamma = Vector3::crossProd(alfa, beta);

	return Matrix33(alfa, beta, gamma);
}


void advanceConstraints(vector<RigidBody>& bodies, float step, bool collisions, float floor){
	int nbodies = bodies.size();

	//Initialize forces with gravity and damping
	vector<Vector3> F(nbodies), T(nbodies);
	for (int i = 0; i < nbodies; i++){
		F[i] = Vector3(0.0f, -bodies[i].Mass()*g, 0.0f);
		F[i] -= d*bodies[i].Velocity();
		T[i] = -d*bodies[i].Omega();
	}

	//Compute free velocities
	for (int i = 0; i<nbodies; i++){
		bodies[i].SetVelocity(bodies[i].Velocity() + (step / bodies[i].Mass()) * F[i]);
		bodies[i].SetOmega(bodies[i].Omega() + step * bodies[i].InertiaInv() * (T[i]+((-1.0f)*Vector3::crossProd(bodies[i].Omega(), bodies[i].Inertia()*bodies[i].Omega()))));
	}

	//Set up joints
	JointFixedConstraint jointFixed(bodies[0], 0, axis, Vector3::ZERO);
	vector<JointConstraint*> joints((int)bodies.size() - 1);
	for (int i = 0; i<nbodies - 1; i++){
		joints[i] = new JointConstraint(bodies[i], bodies[i + 1], i, i + 1, -axis, axis);
	}

	//Check collisions
	vector<ContactConstraint*> contacts;
	if (collisions){
		for (int i = 0; i<nbodies; i++){
			Vector3 positionA = bodies[i].posLocalToGlobal(axis);
			if (positionA[1] < floor){
				contacts.push_back(new ContactConstraint(bodies[i], i, axis, floor));
			}
			Vector3 positionB = bodies[i].posLocalToGlobal(-axis);
			if (positionB[1] < floor){
				contacts.push_back(new ContactConstraint(bodies[i], i, -axis, floor));
			}
		}
	}

	//Set up constraint system
	int cSize = 3 + 3 * joints.size() + contacts.size();
	MatrixMN J(cSize, 6 * nbodies);
	Vector b(cSize);

	//Add the JointCoinstraint J and b for first body. The fixed body
	J.AddBlock33(0, 6 * jointFixed.GetBodyId(), jointFixed.GetJv());
	J.AddBlock33(0, 6 * jointFixed.GetBodyId() + 3, jointFixed.GetJw());
	b.AddBlock3(0, (-1.0f / step)*jointFixed.GetC() - jointFixed.GetV());
	
	//Add J and b for the rest of bodys. JointConstraints
	for (int i = 0; i<joints.size(); i++){
		J.AddBlock33(3 * (i + 1), 6 * joints[i]->GetBodyId1(), joints[i]->GetJv1());
		J.AddBlock33(3 * (i + 1), 6 * joints[i]->GetBodyId1() + 3, joints[i]->GetJw1());
		J.AddBlock33(3 * (i + 1), 6 * joints[i]->GetBodyId2(), joints[i]->GetJv2());
		J.AddBlock33(3 * (i + 1), 6 * joints[i]->GetBodyId2() + 3, joints[i]->GetJw2());
		b.AddBlock3(3 * (i + 1), (-1.0f / step)*joints[i]->GetC() - joints[i]->GetV());
	}

	//Add J and b for the collisionConstraint
	for (int i = 0; i<contacts.size(); i++){
		J.AddBlock13(3 * nbodies + i, 6 * contacts[i]->GetBodyId(), contacts[i]->GetJv());
		J.AddBlock13(3 * nbodies + i, 6 * contacts[i]->GetBodyId() + 3, contacts[i]->GetJw());
		b(3 * nbodies + i) += (-1.0f / step)*contacts[i]->GetC() - contacts[i]->GetV();
	}

	//Set up inverse of mass matrix
	MatrixMN Minv(6 * nbodies, 6 * nbodies);
	for (int i = 0; i<nbodies; i++){
		Minv.AddBlock33(6 * i, 6 * i, (1.0f / bodies[i].Mass())*Matrix33::IDENTITY);
		Minv.AddBlock33(6 * i + 3, 6 * i + 3, bodies[i].InertiaInv());
	}

	//Set up the matrix A for the system
	// A = J * M^-1 * J^T
	MatrixMN A(cSize, cSize);
	A = J * Minv * J.getTranspose();

	//Solve lagrange multipliers
	Vector lambda(cSize);
	if (contacts.size()>0){
		//If collision, resolve LCP using Projected Gauss-Seidel
		bool *flags = new bool[cSize];
		for (int i = 0; i<3 * (joints.size() + 1); i++){
			flags[i] = false;
		}
		for (int i = 3 * (joints.size() + 1); i<cSize; i++){
			flags[i] = true;
		}
		ProjectedGaussSeidel(A, b, lambda, flags);
		delete flags;
	}
	else{
		//If not collision => linear system
		GaussElimination(A, b, lambda);
	}

	//Compute constraint forces loop on constraints and add forces to the bodies
	for (int i = 0; i<nbodies; i++){
		F[i] = Vector3::ZERO;
		T[i] = Vector3::ZERO;
	}
	F[jointFixed.GetBodyId()] += jointFixed.GetJv().getTranspose() * lambda.GetBlock3(0);
	T[jointFixed.GetBodyId()] += jointFixed.GetJw().getTranspose() * lambda.GetBlock3(0);
	for (int i = 0; i<(int)joints.size(); i++){
		F[joints[i]->GetBodyId1()] += joints[i]->GetJv1().getTranspose() * lambda.GetBlock3(3 * (i + 1));
		T[joints[i]->GetBodyId1()] += joints[i]->GetJw1().getTranspose() * lambda.GetBlock3(3 * (i + 1));
		F[joints[i]->GetBodyId2()] += joints[i]->GetJv2().getTranspose() * lambda.GetBlock3(3 * (i + 1));
		T[joints[i]->GetBodyId2()] += joints[i]->GetJw2().getTranspose() * lambda.GetBlock3(3 * (i + 1));
	}
	for (int i = 0; i<(int)contacts.size(); i++){
		F[contacts[i]->GetBodyId()] += contacts[i]->GetJv() * lambda(3 * nbodies + i);
		T[contacts[i]->GetBodyId()] += contacts[i]->GetJw() * lambda(3 * nbodies + i);
	}

	//Compute constrained velocities
	for (int i = 0; i<nbodies; i++){
		bodies[i].SetVelocity(bodies[i].Velocity() + (1.0f / bodies[i].Mass()) * F[i]);
		bodies[i].SetOmega(bodies[i].Omega() + bodies[i].InertiaInv() * T[i]);
	}

	//Integrate positions
	for (int i = 0; i<nbodies; i++){
		bodies[i].SetPosition(bodies[i].Position() + step * bodies[i].Velocity());
		bodies[i].SetRotation(fixRotation(bodies[i].Rotation() + step * Matrix33::ToCrossProd(bodies[i].Omega()) * bodies[i].Rotation()));
		bodies[i].UpdateInertia();
	}

	//Free memory
	for (int i = 0; i<joints.size(); i++){	delete joints[i]; }
	for (int i = 0; i<contacts.size(); i++){ delete contacts[i]; }

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




