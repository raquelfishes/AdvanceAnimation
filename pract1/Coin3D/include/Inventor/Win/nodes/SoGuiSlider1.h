#ifndef COIN_SOGUISLIDER1_H
#define COIN_SOGUISLIDER1_H

// src\Inventor\Win\nodes\SoGuiSlider1.h.  Generated from Slider1.h.in by configure.

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

#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFEnum.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFVec3f.h>

#include <Inventor/nodekits/SoBaseKit.h>

class SoGuiSlider1 : public SoBaseKit {
  typedef SoBaseKit inherited;

  SO_KIT_HEADER(SoGuiSlider1);

  SO_KIT_CATALOG_ENTRY_HEADER(topSeparator);
  SO_KIT_CATALOG_ENTRY_HEADER(surfaceGeometry);
  SO_KIT_CATALOG_ENTRY_HEADER(surfaceTexture);
  SO_KIT_CATALOG_ENTRY_HEADER(surfaceMaterial);
  SO_KIT_CATALOG_ENTRY_HEADER(surfaceTexCoords);
  SO_KIT_CATALOG_ENTRY_HEADER(surfaceCoords);
  SO_KIT_CATALOG_ENTRY_HEADER(surfaceFaceSet);
  SO_KIT_CATALOG_ENTRY_HEADER(knobGeometry);
  SO_KIT_CATALOG_ENTRY_HEADER(knobMaterial);
  SO_KIT_CATALOG_ENTRY_HEADER(knobCoords);
  SO_KIT_CATALOG_ENTRY_HEADER(knobFaceSet);
  SO_KIT_CATALOG_ENTRY_HEADER(knobShadeMaterial);
  SO_KIT_CATALOG_ENTRY_HEADER(knobShadeLineSet);
  SO_KIT_CATALOG_ENTRY_HEADER(knobLightMaterial);
  SO_KIT_CATALOG_ENTRY_HEADER(knobLightLineSet);

public:
  static void initClass(void);

  SoGuiSlider1(void);

  enum Orientation {
    X, Y
  };

  SoSFVec3f size;
  SoSFEnum orientation; // not implemented

  SoSFFloat min;
  SoSFFloat max;
  SoSFFloat value;

  SoSFBool alwaysHook;

  virtual void handleEvent(SoHandleEventAction * action);

  virtual void setSurfaceColor(const SbColor & value);
  virtual void setSurfaceColor(const SbColor & minvalue, const SbColor & maxvalue);
  virtual SbColor getValueAsColor(void) const;

  // FIXME: move to something non-public
  void sizeUpdate(void);
  void valueUpdate(void);
  void minUpdate(void);
  void maxUpdate(void);

protected:
  virtual ~SoGuiSlider1(void);

private:
  void * internals;

};

#endif // !COIN_SOGUISLIDER1_H
