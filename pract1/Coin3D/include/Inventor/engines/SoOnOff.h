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

#ifndef COIN_SOONOFF_H
#define COIN_SOONOFF_H

#include <Inventor/engines/SoSubEngine.h>
#include <Inventor/engines/SoEngineOutput.h>
#include <Inventor/fields/SoSFTrigger.h>


class COIN_DLL_API SoOnOff : public SoEngine {
  typedef SoEngine inherited;

  SO_ENGINE_HEADER(SoOnOff);

public:
  static void initClass(void);
  SoOnOff(void);

  SoSFTrigger on;
  SoSFTrigger off;
  SoSFTrigger toggle;

  SoEngineOutput isOn; // SoSFBool
  SoEngineOutput isOff; // SoSFBool

protected:
  virtual ~SoOnOff(void);

private:
  virtual void evaluate(void);
  virtual void inputChanged(SoField * which);

  SbBool state;
};

#endif // !COIN_SOONOFF_H
