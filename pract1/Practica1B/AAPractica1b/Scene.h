#pragma once


#include <vector>
using namespace std;

#include "Vec3.h"

#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/projectors/SbPlaneProjector.h>


//Forward declarations
class Point;
class Spring;
class GrabForce;
class CGSolver;

class Scene
{
public:
   struct MatrixElem
   {
      float val;
      int row;
      int col;
   };

public:
   // General statics
   static int xPoints;
   static int yPoints;
   static int zPoints;
   static float xSize;
   static float ySize;
   static float zSize;

   static float step;
   static float mass;
   static float stiffness;

protected:
   //Data members
   vector<Point*> points;
   vector<Spring*> springs;

   //Rendering members
   vector<SoSphere*> spheres;
   vector<SoCylinder*> cylinders;
   
   //Transformations
   vector<SoTranslation*> pointXForms;
   vector<SoTransform*> springXForms;

   //Data size
   int nSprings;
   int nPoints;
   int nMovSprings;
   int nMovPoints;

   //Animation
   bool pause;
   bool implicit;
   GrabForce *grabForce;
   SbPlaneProjector *projector;

   //Animation state
   float *x0, *x, *v0, *v;
   Vec3 *v0Block, *vBlock;

   //Conjugate gradient solver
   CGSolver *solver;

public:
   Scene(void);
   Scene(int argc, char* argv[]);
   ~Scene(void);

   //Initialization
   void Init(void);
   void PrintSettings(void);
   void InitAnimation(void);

   //Update scene
   void Update(void);
   void UpdateSprings(void);
   void Animate(void);
   void ComputeForces(void);
   void GetState(float *x, float *v);
   void SetState(float *x, float *v);
   void Integrate(float *x, float *v, float t, float *y);
   void Pause(void);
   void ChangeImplicit(void);
   void Grab(void *data, SbPlaneProjector *proj);
   void UnGrab(void);
   bool Grabbed(void);
   void UpdateGrab(SbVec2f& screenPos);

   //Initialization
   void CreatePoints(void);
   void CreateSprings(void);
   void Add(SoSeparator *root);

};
