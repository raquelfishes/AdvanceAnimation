#include "Scene.h"
#include "Spring.h"
#include "Point.h"
#include "GrabForce.h"
#include "CGSolver.h"

#include <iostream>
using namespace std;

#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSeparator.h>

int Scene::xPoints=5;
int Scene::yPoints=5;
int Scene::zPoints=5;
float Scene::xSize=1.0f;
float Scene::ySize=1.0f;
float Scene::zSize=1.0f;

float Scene::step=0.01f;
float Scene::mass=1.0f;
float Scene::stiffness=300.0f;

Scene::Scene(void)
{
   Init();
   PrintSettings();
}

Scene::Scene(int argc, char* argv[])
{
   int arg=1;
   while(arg<argc)
   {
      //Object size
      if(!strcmp(argv[arg], "-size"))
      {
         xPoints=atoi(argv[++arg]);
         yPoints=atoi(argv[++arg]);
         zPoints=atoi(argv[++arg]);
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
      //Others
      else
      {
         //Print instructions
         cerr << endl << "Unrecognized option " << argv[arg] << endl;
         cerr << "Usage: Test.exe -[option1] [settings] -[option2] [settings] ..." << endl;
         cerr << "Options:" << endl;
         cerr << "\t-method [LeapFrog,BackEuler,RungeKuttaIV]" << endl;
         cerr << "\t-size [sizex] [sizey] [sizez]" << endl;
         cerr << "\t-step [step size in secs]" << endl;
         cerr << "\t-stiff [stiffness value]" << endl;
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
   cerr << "\t-size " << xPoints << " " << yPoints << " " << zPoints << endl;
   cerr << "\t-step " << step << endl;
   cerr << "\t-stiff " << stiffness << endl;
}

void Scene::Init(void)
{
   //Animation settings
   pause=true;
   implicit=true;

   //Create points
   CreatePoints();

   //Create springs
   CreateSprings();

   //Grab force and projector
   projector=NULL;
   grabForce=NULL;

   //Initialize animation
   InitAnimation();
}

void Scene::Pause(void)
{
   pause=!pause;
}

void Scene::ChangeImplicit(void)
{
   implicit=!implicit;
}

void Scene::Update(void)
{
   if(pause)
   {
      return;
   }

   //Perform animation
   Animate();

   //Update the rendering objects
   for(int i=0; i<nPoints; i++)
   {
      //Update sphere xform based on point
      pointXForms[i]->translation.setValue(points[i]->GetPositionRender());
   }
   for(int i=0; i<nSprings; i++)
   {
      //Update cylinder xform based on spring
      springXForms[i]->pointAt(springs[i]->GetPosition(), springs[i]->GetPointA()->GetPositionRender());

      //Update cylinder length
      cylinders[i]->height=springs[i]->GetLength();
   }
   
}

void Scene::Animate(void)
{
   //Get initial state
   GetState(x0, v0);

   //Compute Forces
   ComputeForces();

   //Reset linear system
   solver->ResetSystem();

   //Collect A matrix and b vector from points
   for(int i=0; i<nPoints; i++)
   {
      if(!points[i]->IsFixed())
      {
         points[i]->AddToLinearSystem(solver, implicit);
      }
   }

   //Collect A matrix from springs
   for(int i=0; i<nSprings; i++)
   {
      if(!springs[i]->FullFixed())
      {
         springs[i]->AddToLinearSystem(solver, implicit);
      }
   }

   //Collect A matrix from grabbing force
   if(grabForce!=NULL)
   {
      grabForce->AddToLinearSystem(solver);
   }

   //Solve linear system
   for(int i=0, j=0; i<nMovPoints; i++, j+=3)
   {
      v0Block[i] = Vec3(v0[j], v0[j+1], v0[j+2]);
   }
   solver->SolveCG(vBlock, v0Block);
   for(int i=0, j=0; i<nMovPoints; i++, j+=3)
   {
      v[j] = vBlock[i][0]; v[j+1] = vBlock[i][1]; v[j+2] = vBlock[i][2];
   }

   //Compute positions and set scene state
   Integrate(x0, v, step, x);
   SetState(x, v);
   UpdateSprings();
}

void Scene::UpdateSprings(void)
{
   //Update springs
   for(int i=0; i<nSprings; i++)
   {
      springs[i]->Update();
   }
}

void Scene::GetState(float *x, float *v)
{
   for(int i=0, j=0; i<nPoints; i++)
   {
      if(!points[i]->IsFixed())
      {
         points[i]->GetState(x+j, v+j);
         j+=3;
      }
   }
}

void Scene::SetState(float *x, float *v)
{
   for(int i=0, j=0; i<nPoints; i++)
   {
      if(!points[i]->IsFixed())
      {
         points[i]->SetState(x+j, v+j);
         j+=3;
      }
   }
}

void Scene::Integrate(float *x, float *v, float t, float *y)
{
   for(int i=0; i<nMovPoints*3; i++)
   {
      y[i]=x[i]+t*v[i];
   }
}

void Scene::ComputeForces(void)
{
   //Clear forces
   for(int i=0; i<nPoints; i++)
   {
      points[i]->ClearForce();
   }

   //Add gravity forces
   for(int i=0; i<nPoints; i++)
   {
      points[i]->AddGravity();
   }

   //Add grab force
   if(grabForce!=NULL)
   {
      grabForce->ApplyForce();
   }

   //Add spring forces
   for(int i=0; i<nSprings; i++)
   {
      if(!springs[i]->FullFixed())
      {
         springs[i]->AddForces(implicit);
      }
   }
}

void Scene::CreatePoints(void)
{
   int i, j, k, p;

   //Dimensions
   xSize=1.0f;
   ySize=xSize/xPoints*yPoints;
   zSize=xSize/xPoints*zPoints;

   //Determine number of points
   //x*y*z points on cube corners
   //(x-1)*(y-1)*(z-1) points on cube centers
   int nCorners=xPoints*yPoints*zPoints;
   int nCenters=(xPoints-1)*(yPoints-1)*(zPoints-1);
   nPoints=nCorners+nCenters;

   //Allocate vector
   points.reserve(nPoints);
   for(i=0; i<nPoints; i++)
   {
      points.push_back(new Point());
   }

   //Set point mass
   //Distribute the mass uniformly among all points
   Point::mass=mass/(float)nPoints;

   //Initialize point data
   for(i=0, p=0; i<xPoints; i++)
   {
      for(j=0; j<yPoints; j++)
      {
         for(k=0; k<zPoints; k++, p++)
         {
            points[p]->InitPositionCorner(i, j, k, xPoints, yPoints, zPoints);
         }
      }
   }
   for(i=0; i<xPoints-1; i++)
   {
      for(j=0; j<yPoints-1; j++)
      {
         for(k=0; k<zPoints-1; k++, p++)
         {
            points[p]->InitPositionCenter(i, j, k, xPoints, yPoints, zPoints);
         }
      }
   }

   //Set fixed points and simulation indices
   for(i=0, nMovPoints=0; i<nPoints; i++)
   {
      if(points[i]->GetPosition()[1]>0.9f*ySize)
      {
         points[i]->SetFixed(true);
         points[i]->SetSimIndex(-1);
      }
      else
      {
         points[i]->SetSimIndex(nMovPoints++);
      }
   }

   //Create spheres and transformations
   spheres.reserve(nPoints);
   pointXForms.reserve(nPoints);
   for(i=0; i<nPoints; i++)
   {
      spheres.push_back(new SoSphere());
      spheres[i]->radius.setValue(0.15f/(float)xPoints);
      spheres[i]->setUserData(points[i]);

      pointXForms.push_back(new SoTranslation());
      pointXForms[i]->translation.setValue(points[i]->GetPositionRender());
   }

}

void Scene::CreateSprings(void)
{
   int i, j, k, s, p, q;

   //Determine number of springs
   int nCorners=xPoints*yPoints*zPoints;
   int nXAligned=(xPoints-1)*yPoints*zPoints;
   int nYAligned=xPoints*(yPoints-1)*zPoints;
   int nZAligned=xPoints*yPoints*(zPoints-1);
   int nXYAligned=(xPoints-1)*(yPoints-1)*zPoints;
   int nYZAligned=xPoints*(yPoints-1)*(zPoints-1);
   int nZXAligned=(xPoints-1)*yPoints*(zPoints-1);
   int nCenters=8*(xPoints-1)*(yPoints-1)*(zPoints-1);
   nSprings=nXAligned+nYAligned+nZAligned+nXYAligned+nYZAligned+nZXAligned+nCenters;

   //Allocate vector
   springs.reserve(nSprings);
   for(i=0; i<nSprings; i++)
   {
      springs.push_back(new Spring());
   }

   //Set spring stiffness
   //Distribute the stiffness by unit area in the plane normal to the gravity direction
   Spring::stiffness=stiffness/(float)(xPoints*zPoints);

   //Set spring connections X aligned
   int delta=yPoints*zPoints;
   for(i=0, s=0, p=0; i<xPoints-1; i++)
   {
      for(j=0; j<yPoints; j++)
      {
         for(k=0; k<zPoints; k++, s++, p++)
         {
            springs[s]->Initialize(points[p], points[p+delta]);
         }
      }
   }

   //Set spring connections Y aligned
   delta=zPoints;
   for(i=0, p=0; i<xPoints; i++)
   {
      for(j=0; j<yPoints-1; j++)
      {
         for(k=0; k<zPoints; k++, s++, p++)
         {
            springs[s]->Initialize(points[p], points[p+delta]);
         }
      }
      p+=zPoints;
   }

   //Set spring connections Z aligned
   delta=1;
   for(i=0, p=0; i<xPoints; i++)
   {
      for(j=0; j<yPoints; j++)
      {
         for(k=0; k<zPoints-1; k++, s++, p++)
         {
            springs[s]->Initialize(points[p], points[p+delta]);
         }
         p+=1;
      }
   }

   //Set spring connections XY aligned
   delta=yPoints*zPoints+zPoints;
   for(i=0, p=0; i<xPoints-1; i++)
   {
      for(j=0; j<yPoints-1; j++)
      {
         for(k=0; k<zPoints; k++, s++, p++)
         {
            springs[s]->Initialize(points[p], points[p+delta]);
         }
      }
      p+=zPoints;
   }

   //Set spring connections YZ aligned
   delta=zPoints+1;
   for(i=0, p=0; i<xPoints; i++)
   {
      for(j=0; j<yPoints-1; j++)
      {
         for(k=0; k<zPoints-1; k++, s++, p++)
         {
            springs[s]->Initialize(points[p], points[p+delta]);
         }
         p+=1;
      }
      p+=zPoints;
   }

   //Set spring connections ZX aligned
   delta=yPoints*zPoints+1;
   for(i=0, p=0; i<xPoints-1; i++)
   {
      for(j=0; j<yPoints; j++)
      {
         for(k=0; k<zPoints-1; k++, s++, p++)
         {
            springs[s]->Initialize(points[p], points[p+delta]);
         }
         p+=1;
      }
   }

   //Set spring connections centers
   for(i=0, p=0, q=0; i<xPoints-1; i++)
   {
      for(j=0; j<yPoints-1; j++)
      {
         for(k=0; k<zPoints-1; k++, p++, q++)
         {
            springs[s++]->Initialize(points[nCorners+q], points[p]);
            springs[s++]->Initialize(points[nCorners+q], points[p+yPoints*zPoints]);
            springs[s++]->Initialize(points[nCorners+q], points[p+zPoints]);
            springs[s++]->Initialize(points[nCorners+q], points[p+yPoints*zPoints+zPoints]);
            springs[s++]->Initialize(points[nCorners+q], points[p+1]);
            springs[s++]->Initialize(points[nCorners+q], points[p+yPoints*zPoints+1]);
            springs[s++]->Initialize(points[nCorners+q], points[p+zPoints+1]);
            springs[s++]->Initialize(points[nCorners+q], points[p+yPoints*zPoints+zPoints+1]);
         }
         p+=1;
      }
      p+=zPoints;
   }

   //Set fixed springs
   for(i=0, nMovSprings=0; i<nSprings; i++)
   {
      if(springs[i]->HalfFixed())
      {
         springs[i]->SetSimIndex(-1);
      }
      else
      {
         springs[i]->SetSimIndex(nMovSprings++);
      }
   }

   //Create cylinders and transformations
   cylinders.reserve(nSprings);
   springXForms.reserve(nSprings);
   for(i=0; i<nSprings; i++)
   {
      cylinders.push_back(new SoCylinder());
      cylinders[i]->radius.setValue(0.05f/(float)xPoints);
      cylinders[i]->height.setValue(springs[i]->GetLength());
      cylinders[i]->setUserData(springs[i]);

      springXForms.push_back(new SoTransform());
      springXForms[i]->pointAt(springs[i]->GetPosition(), springs[i]->GetPointA()->GetPositionRender());
   }

}

void Scene::InitAnimation(void)
{
   x=new float[3*nMovPoints];
   v=new float[3*nMovPoints];
   x0=new float[3*nMovPoints];
   v0=new float[3*nMovPoints];

   vBlock=new Vec3[nMovPoints];
   v0Block=new Vec3[nMovPoints];

   //Initialize solver
   solver = new CGSolver(nMovPoints, nMovSprings);
   solver->Init();
   solver->SetMaxIter(10);
   solver->SetTolerance(1e-10f);

   //Set indices of sparse blocks
   for (int i=0; i<nSprings; i++)
   {
      springs[i]->InitLinearSystem(solver);
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
      Point *aux=points[i];
      if(!points[i]->IsFixed())
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
      if(points[i]->IsFixed())
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

void Scene::Grab(void *data, SbPlaneProjector *proj)
{
   projector=proj;
   grabForce=new GrabForce((Point*)data);
}

void Scene::UnGrab(void)
{
   delete projector;
   projector=NULL;

   delete grabForce;
   grabForce=NULL;
}

bool Scene::Grabbed(void)
{
   return grabForce!=NULL;
}

void Scene::UpdateGrab(SbVec2f& screenPos)
{
   SbVec3f pos = projector->project(screenPos);
   grabForce->UpdateGrab(Vec3(pos[0], pos[1], pos[2]));
}

