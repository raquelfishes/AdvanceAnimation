#ifndef COIN_SOTEXTUREQUALITYELEMENT_H
#define COIN_SOTEXTUREQUALITYELEMENT_H

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

#include <Inventor/elements/SoFloatElement.h>

class COIN_DLL_API SoTextureQualityElement : public SoFloatElement {
  typedef SoFloatElement inherited;

  SO_ELEMENT_HEADER(SoTextureQualityElement);
public:
  static void initClass(void);
protected:
  virtual ~SoTextureQualityElement();

public:
  virtual void init(SoState * state);
  static void set(SoState * const state, SoNode * const node,
                  const float quality);
  static void set(SoState * const state, const float quality);
  static float get(SoState * const state);
  static float getDefault(void);
};

#ifndef COIN_INTERNAL
#include <Inventor/elements/SoTextureOverrideElement.h>
#endif // !COIN_INTERNAL

#endif // !COIN_SOTEXTUREQUALITYELEMENT_H
