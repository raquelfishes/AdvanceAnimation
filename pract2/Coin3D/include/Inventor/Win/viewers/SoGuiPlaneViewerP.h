#ifndef SOGUIPLANEVIEWERP_H
#define SOGUIPLANEVIEWERP_H

// src\Inventor\Win\viewers\SoGuiPlaneViewerP.h.  Generated from SoGuiPlaneViewerP.h.in by configure.

/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) by Kongsberg Oil & Gas Technologies.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Kongsberg Oil & Gas Technologies
 *  about acquiring a Coin Professional Edition License.
 *
 *  See http://www.coin3d.org/ for more information.
 *
 *  Kongsberg Oil & Gas Technologies, Bygdoy Alle 5, 0257 Oslo, NORWAY.
 *  http://www.sim.no/  sales@sim.no  coin-support@coin3d.org
 *
\**************************************************************************/

#ifndef SOWIN_INTERNAL
#error this is a private header file
#endif /* !SOWIN_INTERNAL */

#include <Inventor/SbLinear.h>

class SoWinPlaneViewer;

// ************************************************************************

// This class contains private data and methods used within the
// SoGuiPlaneViewer class.

class SoGuiPlaneViewerP
{
public:
  ~SoGuiPlaneViewerP();

  void commonConstructor(void);

  void pan(const SbVec2f & thispos, const SbVec2f & prevpos);
  void rotateZ(const float angle) const;

  void viewPlaneX(void) const;
  void viewPlaneY(void) const;
  void viewPlaneZ(void) const;

  void setCanvasSize(const SbVec2s size);
  void setPointerLocation(const SbVec2s location);
  int getPointerXMotion(void) const;
  int getPointerYMotion(void) const;
  float getPointerOrigoAngle(void) const;
  float getPointerOrigoMotionAngle(void) const;

  void updateAnchorScenegraph(void) const;

  enum PlaneViewerMode {
    SCENEGRAPH_INTERACT_MODE,

    IDLE_MODE,

    DOLLY_MODE,
    TRANSLATE_MODE,

    ROTZ_WAIT_MODE,
    ROTZ_MODE,

    SEEK_WAIT_MODE,
    SEEK_MODE
  } mode;

  void changeMode(PlaneViewerMode newmode);
  void setCursorRepresentation(PlaneViewerMode mode);

  struct pointerdata {
    SbVec2s now;
    SbVec2s then;
  } pointer;
  SbVec2s canvas;

  SbBool ctrldown;
  SbBool shiftdown;
  SbBool button1down;
  SbBool button3down;

  SbPlane panningplane;

  class SoDrawStyle * lineds[2];

  class SoNode * superimposition;
  struct superdata {
    class SoCoordinate3 * coords;
    class SoOrthographicCamera * camera;
  } super;

protected:
  SoGuiPlaneViewerP(SoWinPlaneViewer * publ);
  SoWinPlaneViewer * pub;
};

// ************************************************************************

#endif // ! SOGUIPLANEVIEWERP_H
