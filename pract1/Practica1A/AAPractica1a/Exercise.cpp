#include <iostream>
using namespace std;
#include "Vec2.h"
#include "LinSys.h"
#include "Exercise.h"

// gravitational acceleration (9.81)
const float g = 9.81f;

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

class Spring1D
{
public:
   float p1;
   float p2;
   float L;
   float k;

   Spring1D(void){}
   Spring1D(float p1, float p2, float L, float k) : p1(p1), p2(p2), L(L), k(k) {}
   ~Spring1D(void){}

   float force1(void) const { return -k*(p1-p2-L); }
   float force2(void) const { return -k*(p2-p1+L); }

   float dF1dp1(void) const { return -k; }
   float dF1dp2(void) const { return k; }
   float dF2dp1(void) const { return k; }
   float dF2dp2(void) const { return -k; }

};

class Spring2D
{
public:
   Vec2 p1;
   Vec2 p2;
   float L0;
   float L;
   float k;
   Vec2 u1;
   Vec2 F;
   Matrix2 dFdp;

   Spring2D(void){}
   Spring2D(const Vec2& p1, const Vec2& p2, float L0, float k) : p1(p1), p2(p2), L0(L0), k(k)
   {
      L = (p1-p2).length();
      u1 = (1.0f/L) * (p1 - p2);
      F = (k*(L0-L)) * u1;

      //u1 = (p1 - p2)/l
      //l^2 = (p1 - p2)^T (p1 - p2)
      //2 l dl/dp1 = 2 (p1 - p2)^T
      //dl/dp1 = u1^T

      //l * u1 = (p1 - p2)
      //l * du1/dp1 + u1 * dl/dp1 = I
      //du1/dp1 = 1/l * (I - u1 * dl/dp1)
      //du1/dp1 = 1/l * (I - u1 * u1^T)

      //dF1/dp1 = - k * (l - l0) * du1/dp1 - k * u1 * dl/dp1
      //dF1/dp1 = - k * (l - l0) * 1/l * (I - u1 * u1^T) - k * u1 * u1^T
      //dF1/dp1 = - k * (1 - l0/l) * (I - u1 * u1^T) - k * u1 * u1^T
      //dF1/dp1 = - k * (I - u1 * u1^T) + k * l0/l * (I - u1 * u1^T) - k * u1 * u1^T
      //dF1/dp1 = - k * I + k * l0/l * (I - u1 * u1^T)
      //dF1/dp1 = - k * (1 - l0/l) * I - k * l0/l * u1 * u1^T

      dFdp = (k*(L0/L-1)) * Matrix2::IDENTITY - ((k*L0/L) * u1) * u1;
   }
   ~Spring2D(void){}

   const Vec2& force1(void) const { return F; }
   Vec2 force2(void) const { return (-1.0)*F; }

   const Matrix2& dF1dp1(void) const { return dFdp; }
   Matrix2 dF1dp2(void) const { return (-1.0)*dFdp; }
   Matrix2 dF2dp1(void) const { return (-1.0)*dFdp; }
   const Matrix2& dF2dp2(void) const { return dFdp; }

};


class Triangle
{
public:
   //3 points defined in counterclockwise order: p1, p2, p3
   Vec2 p1;
   Vec2 p2;
   Vec2 p3;
   float A;
   float A0;
   float k;
   Vec2 dAdp1T, dAdp2T, dAdp3T;
   Matrix2 d2AdpTdnext, d2AdpTdprev;

   Triangle(void){}
   Triangle(const Vec2& p1, const Vec2& p2, const Vec2& p3, float A0, float k) : p1(p1), p2(p2), p3(p3), A0(A0), k(k)
   {
      //A = 1/2 * (p2 - p1) x (p3 - p1) = 1/2 * (p3 - p2) x (p1 - p2) = 1/2 * (p1 - p3) x (p2 - p3)
      A = 0.5f * (p2 - p1).cross(p3 - p1);

      //A = 1/2 * (prev - next) x (p - next)
      //dAdp = 1/2 * (prev - next)*
      dAdp1T = 0.5f * Vec2::CrossProd(p3 - p2);
      dAdp2T = 0.5f * Vec2::CrossProd(p1 - p3);
      dAdp3T = 0.5f * Vec2::CrossProd(p2 - p1);

      //u* = (-y, x)
      //du*T/du = (0, -1; 1, 0)
      d2AdpTdprev = Matrix2(0.0f, -1.0f, 1.0f, 0.0f);
      d2AdpTdnext = Matrix2(0.0f, 1.0f, -1.0f, 0.0f);
      //d2AdpTdp = 0

      //E = 1/2 * k * (A - A0)^2

      //Fi = - k * (A - A0) * dA/dpi^T

      //dFi/dpj = - k * (A - A0) * d(dA/dpi^T)/dpj - k * dA/dpi^T * dA/dpj
   }
   ~Triangle(void){}

   Vec2 force1(void) const { return (k*(A0-A)) * dAdp1T; }
   Vec2 force2(void) const { return (k*(A0-A)) * dAdp2T; }
   Vec2 force3(void) const { return (k*(A0-A)) * dAdp3T; }

   Matrix2 dF1dp1(void) const { return ((-k) * dAdp1T) * dAdp1T; }
   Matrix2 dF1dp2(void) const { return (k*(A0-A)) * d2AdpTdnext - (k * dAdp1T) * dAdp2T; }
   Matrix2 dF1dp3(void) const { return (k*(A0-A)) * d2AdpTdprev - (k * dAdp1T) * dAdp3T; }
   Matrix2 dF2dp1(void) const { return (k*(A0-A)) * d2AdpTdprev - (k * dAdp2T) * dAdp1T; }
   Matrix2 dF2dp2(void) const { return ((-k) * dAdp2T) * dAdp2T; }
   Matrix2 dF2dp3(void) const { return (k*(A0-A)) * d2AdpTdnext - (k * dAdp2T) * dAdp3T; }
   Matrix2 dF3dp1(void) const { return (k*(A0-A)) * d2AdpTdnext - (k * dAdp3T) * dAdp1T; }
   Matrix2 dF3dp2(void) const { return (k*(A0-A)) * d2AdpTdprev - (k * dAdp3T) * dAdp2T; }
   Matrix2 dF3dp3(void) const { return ((-k) * dAdp3T) * dAdp3T; }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Exercise 1
// hanging mass point
void AdvanceTimeStep1(float k, float m, float d, float L, float kc, float dt, int meth, float p1, float v1, float& p2, float& v2, float& p2old, bool collision)
{
	switch(meth){
		case 1:
			// Euler Explicit
			eulerExplicit(k,m,d,L,kc,dt,p1,v1,p2,v2,collision);
			break;
		case 2:
			// Euler Symplectic
			eulerSymplectic(k,m,d,L,kc,dt,p1,v1,p2,v2,collision);
			break;
		case 3:
			// Mid Point
			midPoint(k,m,d,L,kc,dt,p1,v1,p2,v2,collision);
			break;
		case 4:
			// Euler Implicit
			eulerImplicit(k,m,d,L,kc,dt,p1,v1,p2,v2,collision);
			break;
		case 5:
			// Verlet
			verlet(k,m,d,L,kc,dt,p1,v1,p2,v2,p2old,collision);
			break;
		default:
			break;
	}
}


// Exercise 2
// square
void AdvanceTimeStep2(float k, float m, float d, float L, float kA, float A, float dt, int method,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area)
{
	switch(method){
		case 2:
			// Euler Symplectic
			eulerSymplectic2(k,m,d,L,kA,A,dt,p1,v1,p2,v2,p3,v3,p4,v4,springs,area);
			break;
		case 4:
			// Euler Implicit
			//eulerImplicit2(k,m,d,L,kc,dt,p1,v1,p2,v2,collision);
			break;
		default:
			break;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Euler explicit method
void eulerExplicit(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision){
	//Forces
	float mg = m*(-g);
	float f2 = 0.0f;

	//Add weight force to f2
	f2 += mg;

	//Add Springs forces
	Spring1D sp1(p1,p2,L,k);
	f2 += sp1.force2();
	if (p2<0.0f && collision){
		Spring1D sp2(0.0f,p2,0.0f,kc);
		f2 += sp2.force2();
	}

	//Add dumping forces
	f2-=d*v2;

	//Integrate
	//Position uses v prev
	p2+=dt*v2;
	v2+=dt*(1.0f/m)*f2;

}

void eulerSymplectic(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision){
	//Forces
	float mg = m*(-g);
	float f2 = 0.0f;

	//Add weight force to f2
	f2 += mg;

	//Add Spring forces
	Spring1D sp1(p1,p2,L,k);
	f2 += sp1.force2();
	if (p2<0.0f && collision){
		Spring1D sp2(0.0f,p2,0.0f,kc);
		f2 += sp2.force2();
	}

	//Add dumping forces
	f2-=d*v2;

	//Integrate
	//Position uses v new
	v2+=dt*(1.0f/m)*f2;
	p2+=dt*v2;
}

void midPoint(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision){
	//Forces
	float mg = m*(-g);

	//Calculate aceleration: a0 = F/m at t
	Spring1D sp1(p1,p2,L,k);
	float a0 = (1.0f/m)*(mg - d*v2 + sp1.force2());
	if (p2<0.0f && collision){
		Spring1D sp2(0.0f,p2,0.0f,kc);
		a0-=(1.0f/m)*(sp2.force2());
	}

	//Calculate velocity: v2 at dt/2
	float v2_aux = v2+(dt/2.0f)*a0;

	//Calculate position: p2 at dt/2
	float p2_aux = p2+(dt/2.0f)*v2_aux;

	//Calculate aceleration: a2 at dt/2
	sp1 = Spring1D(p1,p2_aux,L,k);
	float a2 = (1.0f/m)*(mg-d*v2_aux*sp1.force2());
	/****** PREGUNTAR A OTADUY, HAY QUE VOLVER A CALCULAR COLISION??? ******/
	//if (p2_aux<0.0f && collision){
	//	Spring1D sp2(0.0f,p2,abs(p2),kc);
	//	a0-=(1.0f/m)*(sp2.force2());
	//}

	//Integrate
	p2+=dt*v2_aux;
	v2+=dt*a2;
}

void eulerImplicit(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, bool collision){
	//Forces
	float mg = m*(-g);
	float f2 = 0.0f;

	//Add Spring forces
	Spring1D sp1(p1,p2,L,k);
	f2 += mg - d*v2 + sp1.force2();
	if (p2<0.0f && collision){
		Spring1D sp2(0.0f,p2,0.0f,kc);
		f2 += sp2.force2();
	}

	//Set up system
	//v(h)=v(0)+h/m*F(0)+h/m*dFdv*v(0) / 1-h/m*dF/dv-h^2/m*dF/dx
	
	//b=v(0)+h/m*F(0)+h/m*dFdv*v(0) 
	//We don't want have h/m so we multiplicate with m
	//b=v(0)*m+h*F(0)+h*dFdv*v(0)
	float b = (m+dt*d)*v2 + dt*f2;
	
	//m_impl = 1-h/m*dF/dv-h^2/m*dF/dx
	//We don't want have h/m so we multiplicate with m
	//m_impl = m-h*dF/dv-h^2*dF/dx
	float m_impl = (m+dt*d-dt*dt*sp1.dF2dp2());

   //Solve for velocity
   v2 = b/m_impl;

   //Integrate position
   p2+=dt*v2;

}

void verlet(float k, float m, float d, float L, float kc, float dt, float p1, float v1, float& p2, float& v2, float& p2old, bool collision){
	//Saving old position
	float p_aux = p2;

	//Forces
	float mg = m*(-g);
	float f2 = 0.0f;
	
	//Add weight force to f2
	f2 += mg;

	//Add Spring forces
	Spring1D sp1(p1,p2,L,k);
	f2 += sp1.force2();
	if (p2<0.0f && collision){
		Spring1D sp2(0.0f,p2,abs(p2),kc);
		f2 += sp2.force2();
	}

	//Add damping forces
	f2-=d*v2;

	//Integrate
	float a = (1.0f/m)*f2;
	p2=2.0f*p2-p2old+dt*dt*a;
	v2=(p2-p2old)/(2.0f*dt);
	
	//Updating old position
	p2old=p_aux;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void eulerSymplectic2(float k, float m, float d, float L, float kA, float A, float dt,
   const Vec2& p1, const Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3, Vec2& p4, Vec2& v4, bool springs, bool area){
	
	cout << A << endl;

	//Initialize forces
	Vec2 mg(0.0f, -m*g);

	//Vec2 f1(0.0f, 0.0f);
	Vec2 f2(0.0f, 0.0f);
	Vec2 f3(0.0f, 0.0f);
	Vec2 f4(0.0f, 0.0f);

	//Add weight force
	//f1+=mg;
	f2+=mg;
	f3+=mg;
	f4+=mg;

	//Add spring forces
   if (springs)
   {
      Spring2D sp12(p1, p2, L, k);
      Spring2D sp23(p2, p3, L, k);
	  Spring2D sp34(p3, p4, L, k);
      Spring2D sp41(p4, p1, L, k);
      //f1+=sp12.force1()+sp41.force2();
      f2+=sp23.force1()+sp12.force2();
      f3+=sp34.force1()+sp23.force2();
	  f4+=sp41.force1()+sp34.force2();
   }

   //Add damping forces
   //f1-=d*v1;
   f2-=d*v2;
   f3-=d*v3;
   f4-=d*v4;
   
    //Add triangle forces

   if (area)
   {
      Triangle tri1(p1, p2, p3, A/2, kA);
	  Triangle tri2(p3, p4, p1, A/2, kA);
      //f1+=tri1.force1()+tri2.force3();
      f2+=tri1.force2();
      f3+=tri1.force3()+tri2.force1();
	  f4+=tri2.force2();
   }
   
   //Integration
   v2+=dt*(1.0f/m)*f2;
   p2+=dt*v2;
   v3+=dt*(1.0f/m)*f3;
   p3+=dt*v3;
   v4+=dt*(1.0f/m)*f4;
   p4+=dt*v4;

   cout << "p2: " << p2.x << " " << p2.y << " v2: " << v2.x << " " << v2.y << "\np3: " << p3.x << " " << p3.y << " v3: " << v3.x << " " << v3.y << "\np4: " << p4.x << " " << p4.y << " v4: " << v4.x << " " << v4.y << endl;
}


