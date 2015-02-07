// coin3D includes
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>

#include <Inventor/sensors/SoIdleSensor.h> 
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/events/SoLocation2Event.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/projectors/SbPlaneProjector.h>
#include <Inventor/nodes/SoCamera.h>

#include ".\Scene.h"
#include <iostream>


struct StrAllData
{
   SoWinExaminerViewer *viewer;
   Scene *scene;
};
typedef struct StrAllData AllData ;


//Idle function
static void idlefunc(void *data, SoSensor* s)
{
   Scene *scene=((AllData*)data)->scene;

   scene->Update();

   s->schedule();
}

//Keyboard callback
static void keyboardfunc(void *data, SoEventCallback *eventCB)
{
   Scene *scene=((AllData*)data)->scene;
   const SoEvent *event = eventCB->getEvent();
   if(SO_KEY_PRESS_EVENT(event, P))
   {
      scene->Pause();
      eventCB->setHandled();
   }
   else if(SO_KEY_PRESS_EVENT(event, S))
   {
      scene->Springs();
      eventCB->setHandled();
   }
   else if(SO_KEY_PRESS_EVENT(event, A))
   {
      scene->Area();
      eventCB->setHandled();
   }
   else if(SO_KEY_PRESS_EVENT(event, C))
   {
      scene->Collision();
      eventCB->setHandled();
   }
}

int main(int argc, char* argv[])
{
   HWND window = SoWin::init(argv[0]);
   if (window==NULL) exit(1);

   SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);

   //Create root
   SoSeparator * root = new SoSeparator;
   root->ref();

   //Initialize scene and add to root
   Scene *scene=new Scene(argc, argv);
   scene->Add(root);

   //Create application data structure
   AllData appData;
   appData.scene=scene;
   appData.viewer=viewer;

   // Create idle sensor and pass pointer to scene
   SoIdleSensor *idle = new SoIdleSensor(idlefunc, &appData);
   idle->schedule();

   // Create event handler for keyboard
   SoEventCallback *keyboardEvent = new SoEventCallback;
   keyboardEvent->addEventCallback(SoKeyboardEvent::getClassTypeId(), keyboardfunc, &appData);
   root->addChild(keyboardEvent);

   //Start show
   viewer->setSceneGraph(root);
   viewer->show();

   SoWin::show(window);
   SoWin::mainLoop();
   delete viewer;
   root->unref();
   return 0;
}

