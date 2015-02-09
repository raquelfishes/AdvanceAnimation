#include ".\scene.h"
#include ".\primitives.h"
#include ".\exercise.h"
#include "Vec2.h"

#include <iostream>
using namespace std;

#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/SbTime.h>

float Scene::step=0.01f;
float Scene::mass=1.0f;
float Scene::stiffness=10.0f;
float Scene::stiffpenalty=100.0f;
float Scene::stiffarea=10.0f;
float Scene::damping=0.01f;

Scene::Method Scene::method=EULER_SYMPLEC;
Scene::Testcase Scene::testcase=SPRING1D;

Scene::Scene(void)
{
   Init();
   PrintSettings();
}

//some default call: -testcase hanging -method Euler -stiff 10 -mass 0.1 -step 0.003 -damp 0.01

Scene::Scene(int argc, char* argv[])
{
   // defaults:
   //testcase  = SPRING1D;
   ////method    = verlet;
   ////method    = euler_symplec;
   //method    = EULER;
   //stiffness = 1.0f;
   //stiffpenalty = 100.0f;
   //mass      = 0.1f;
   //step      = 0.003f;
   //damping   = 0.01f;
   testcase  = SQUARE;
   //method    = BACK_EULER;
   method    = EULER_SYMPLEC;
   stiffness = 100.0f;
   stiffarea = 100.0f;
   mass      = 0.1f;
   step      = 0.025f;
   damping   = 0.1f;

   int arg=1;
   while(arg<argc)
   {
      //Testcase
      if(!strcmp(argv[arg], "-testcase"))
      {
         arg++;
         //for(unsigned int i=0; i<strlen(argv[arg]); i++) argv[arg][i] = tolower(argv[arg][i]);
         if(!strcmp(argv[arg], "square"))
         {
            testcase=SQUARE;
         }
         else if(!strcmp(argv[arg], "spring1d"))
         {
            testcase=SPRING1D;
         }
         else
         {
            cerr << "Unrecognized testcase " << argv[arg] << endl;
            exit(1);
         }
         arg++;
      }
      //Integration method
      else if(!strcmp(argv[arg], "-method"))
      {
         arg++;
         //for(unsigned int i=0; i<strlen(argv[arg]); i++) argv[arg][i] = tolower(argv[arg][i]);
         if(!strcmp(argv[arg], "euler"))
         {
            method=EULER;
         }
         else if(!strcmp(argv[arg], "backeuler"))
         {
            method=BACK_EULER;
         }
         else if(!strcmp(argv[arg], "eulersymplec"))
         {
            method=EULER_SYMPLEC;
         }
         else if(!strcmp(argv[arg], "midpoint"))
         {
            method=MIDPOINT;
         }
         else if(!strcmp(argv[arg], "verlet"))
         {
            method=VERLET;
         }
         else
         {
            cerr << "Unrecognized method " << argv[arg] << endl;
            exit(1);
         }
         arg++;
      }
      //Step size
      else if(!strcmp(argv[arg], "-step"))
      {
         step=(float)atof(argv[++arg]);
         arg++;
      }
      //Stiffness
      else if(!strcmp(argv[arg], "-stiff"))
      {
         stiffness=(float)atof(argv[++arg]);
         arg++;
      }
      //Penalty
      else if(!strcmp(argv[arg], "-penalty"))
      {
         stiffpenalty=(float)atof(argv[++arg]);
         arg++;
      }
      //Area preservation
      else if(!strcmp(argv[arg], "-area"))
      {
         stiffarea=(float)atof(argv[++arg]);
         arg++;
      }
      //Damping
      else if(!strcmp(argv[arg], "-damp"))
      {
         damping=(float)atof(argv[++arg]);
         arg++;
      }
      //Mass
      else if(!strcmp(argv[arg], "-mass"))
      {
         mass=(float)atof(argv[++arg]);
         arg++;
      }
      //Others
      else
      {
         //Print instructions
         cerr << endl << "Unrecognized option " << argv[arg] << endl;
         cerr << "Usage: Test.exe -[option1] [settings] -[option2] [settings] ..." << endl;
         cerr << "Options:" << endl;
         cerr << "\t-testcase [spring1d,square]" << endl;
         cerr << "\t-method [euler,eulersymplec,midpoint,backeuler,verlet]" << endl;
         cerr << "\t-step [step size in secs]" << endl;
         cerr << "\t-stiff [spring stiffness]" << endl;
         cerr << "\t-penalty [stiffness value for penalty contact (Exercise 2)]" << endl;
         cerr << "\t-area [stiffness value for area preservation (Exercise 2)]" << endl;
         cerr << "\t-damp [damping value]" << endl << endl;
         break;
      }
   }

   PrintSettings();
   Init();
}

Scene::~Scene(void)
{
}

void Scene::PrintSettings(void)
{
   cerr << endl << "Current Settings:" << endl;
   cerr << "\t-testcase " << (int)testcase << endl;
   if(method==BACK_EULER)
   {
      cerr << "\t-method BackEuler" << endl;
   }
   else if(method==MIDPOINT)
   {
      cerr << "\t-method Midpoint" << endl;
   }
   else if(method==EULER)
   {
      cerr << "\t-method Euler" << endl;
   }
   else if(method==EULER_SYMPLEC)
   {
      cerr << "\t-method EulerSymplec" << endl;
   }
   else if(method==VERLET)
   {
      cerr << "\t-method Verlet" << endl;
   }
   else 
   {
      cerr << "\t-method ?" << endl;
   }
   cerr << "\t-mass " << mass << endl;
   cerr << "\t-step " << step << endl;
   cerr << "\t-stiff " << stiffness << endl;
   cerr << "\t-penalty " << stiffpenalty << endl;
   cerr << "\t-area " << stiffarea << endl;
   cerr << "\t-damp " << damping << endl << endl;
}

void Scene::Init(void)
{
   //Animation settings
   pause=true;
   collision=true;
   area=true;
   length=true;

   //Create points & springs
   nPoints = testcase == SPRING1D ? 2 : 4;

   //Allocate vector
   points.reserve(nPoints);
   for(int i=0; i<nPoints; i++)
   {
      points.push_back(new Sphere());
   }

   //Set points at rest configuration
   const float pi=3.1416f;
   p1 = Vec2(0, 1);
   p2old = p2 = Vec2(0, 0.5f);
   p3 = Vec2(0.5f, 0.5f);
   p4 = Vec2(0.5f, 1.0f);
   v1 = v2 = v3 = Vec2::ZERO;
   L = 0.5f;
   A = 2 * 0.5f * 0.5f;

   points[0]->position = SbVec3f(p1.x,p1.y, 0.);
   points[1]->position = SbVec3f(p2.x,p2.y, 0.);
   if(nPoints > 2) {
      points[2]->position = SbVec3f(p3.x,p3.y, 0.);
      points[3]->position = SbVec3f(p4.x,p4.y, 0.);
   }
   points[0]->fixed = true;

   //Create spheres and transformations
   spheres.reserve(nPoints);
   pointXForms.reserve(nPoints);
   for(int i=0; i<nPoints; i++)
   {
      spheres.push_back(new SoSphere());
      spheres[i]->radius.setValue(0.08f);
      spheres[i]->setUserData(points[i]);

      pointXForms.push_back(new SoTranslation());
      pointXForms[i]->translation.setValue(points[i]->position);
   }

   nSprings = nPoints == 2 ? 1 : 4;

   //Allocate vector
   springs.reserve(nSprings);
   for(int i=0; i<nSprings; i++)
   {
      springs.push_back(new Tube());
   }
   springs[0]->Initialize(points[0],points[1]);
   if (nPoints == 4)
   {
      springs[1]->Initialize(points[1],points[2]);
      springs[2]->Initialize(points[2],points[3]);
      springs[3]->Initialize(points[3],points[0]);
   }

   //Create cylinders and transformations
   cylinders.reserve(nSprings);
   springXForms.reserve(nSprings);
   for(int i=0; i<nSprings; i++)
   {
      cylinders.push_back(new SoCylinder());
      cylinders[i]->radius.setValue(0.03f);
      cylinders[i]->height.setValue(springs[i]->length);
      cylinders[i]->setUserData(springs[i]);

      springXForms.push_back(new SoTransform());
      springXForms[i]->pointAt(springs[i]->position, springs[i]->x->position);
   }
}


void Scene::Update(void)
{
   if(pause)
   {
      return;
   }

   //Perform animation
   switch(testcase) {
   case SPRING1D:
      AdvanceTimeStep1(stiffness, mass, damping, L, stiffpenalty, step, method, p1.y, v1.y, p2.y, v2.y, p2old.y, collision);
      break;
   case SQUARE:
      AdvanceTimeStep2(stiffness, mass, damping, L, stiffarea, A, step, method, p1, v1, p2, v2, p3, v3, p4, v4, length, area);
      break;
   }
   points[0]->position = SbVec3f(p1.x,p1.y, 0.);
   points[1]->position = SbVec3f(p2.x,p2.y, 0.);
   if(nPoints>2) {
      points[2]->position = SbVec3f(p3.x,p3.y, 0.); 
      points[3]->position = SbVec3f(p4.x,p4.y, 0.); 
   }
   for(int i=0; i<nSprings; i++) springs[i]->Update();

   //Update the rendering objects
   for(int i=0; i<nPoints; i++)
   {
      //Update sphere xform based on point
      pointXForms[i]->translation.setValue(points[i]->position);
   }
   for(int i=0; i<nSprings; i++)
   {
      //Update cylinder xform based on spring
      springXForms[i]->pointAt(springs[i]->position, springs[i]->x->position);

      //Update cylinder length
      cylinders[i]->height=springs[i]->length;
   }

   // reduce frame rate a bit
   SbTime t = SbTime::getTimeOfDay();
   while(SbTime::getTimeOfDay().getValue()-t.getValue()<0.00125) {
      // delay
   }
}




void Scene::Add(SoSeparator *root)
{
   //Add separator for moving points
   SoSeparator *pointSeparator=new SoSeparator;
   root->addChild(pointSeparator);

   //Moving point properties
   SoMaterial *pointMaterial=new SoMaterial;
   pointMaterial->diffuseColor.setValue(0.9f, 0.0f, 0.0f);
   pointSeparator->addChild(pointMaterial);

   //Add points with xforms
   for(int i=0; i<nPoints; i++)
   {
      Sphere *aux=points[i];
      if(!points[i]->fixed)
      {
         SoSeparator *pSeparator=new SoSeparator;
         pointSeparator->addChild(pSeparator);
         pSeparator->addChild(pointXForms[i]);
         pSeparator->addChild(spheres[i]);
      }
   }

   //Add separator for fixed points
   SoSeparator *pointFixedSeparator=new SoSeparator;
   root->addChild(pointFixedSeparator);

   //Fixed point properties
   SoMaterial *pointMaterialFixed=new SoMaterial;
   pointMaterialFixed->diffuseColor.setValue(0.0f, 0.0f, 0.9f);
   pointFixedSeparator->addChild(pointMaterialFixed);

   //Add points with xforms
   for(int i=0; i<nPoints; i++)
   {
      if(points[i]->fixed)
      {
         SoSeparator *pSeparator=new SoSeparator;
         pointFixedSeparator->addChild(pSeparator);
         pSeparator->addChild(pointXForms[i]);
         pSeparator->addChild(spheres[i]);
      }
   }

   //Add separator for springs
   SoSeparator *springSeparator=new SoSeparator;
   root->addChild(springSeparator);

   //Spring properties
   SoMaterial *springMaterial=new SoMaterial;
   springMaterial->diffuseColor.setValue(0.9f, 0.9f, 0.9f);
   springSeparator->addChild(springMaterial);

   //Create common rotation
   SoTransform *cylOrient=new SoTransform;
   cylOrient->pointAt(SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(0.0f, -1.0f, 0.0f));

   //Add springs with xforms
   for(int i=0; i<nSprings; i++)
   {
      SoSeparator *sSeparator=new SoSeparator;
      springSeparator->addChild(sSeparator);
      sSeparator->addChild(springXForms[i]);
      sSeparator->addChild(cylOrient);
      sSeparator->addChild(cylinders[i]);
   }
}

