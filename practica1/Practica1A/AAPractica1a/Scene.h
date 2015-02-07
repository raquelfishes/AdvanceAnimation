#pragma once

#include <vector>
using namespace std;

#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/projectors/SbPlaneProjector.h>

#include "Vec2.h"

//Forward declarations
class Sphere;
class Tube;

class Scene
{

public:
   // Scene parameters
   static float step;
   static float mass;
   static float stiffness;
   static float stiffpenalty;
   static float stiffarea;
   static float damping;
   float L;
   float A;
   Vec2 p1,p2,p3,p4,p2old;
   Vec2 v1,v2,v3,v4;

   static enum Method{EULER=1, EULER_SYMPLEC=2, MIDPOINT=3, BACK_EULER=4, VERLET=5};
   static Method method;
   static enum Testcase{SPRING1D=1,SQUARE=2};
   static Testcase testcase;

protected:
   //Data members
   vector<Sphere*> points;
   vector<Tube*> springs;

//Rendering members
   vector<SoSphere*> spheres;
   vector<SoCylinder*> cylinders;
   
   //Transformations
   vector<SoTranslation*> pointXForms;
   vector<SoTransform*> springXForms;

   //Data size
   int nSprings;
   int nPoints;

   //Animation
   bool collision;
   bool area;
   bool length;
   bool pause;

   //Animation state
   float *x0, *x;
   float *v0, *v;


public:
   Scene(void);
   Scene(int argc, char* argv[]);
   ~Scene(void);

   //Initialization
   void Init(void);
   void PrintSettings(void);
   void Add(SoSeparator *root);

   //Update scene
   void Update(void);
   void Animate(void);
   void Pause(void) { pause=!pause; }
   void Springs(void) { length=!length; }
   void Area(void) { area=!area; }
   void Collision(void) { collision=!collision; }
};
