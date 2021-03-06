#ifndef COIN_SBPLIST_H
#define COIN_SBPLIST_H

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

#ifndef COIN_INTERNAL
// The next includes are for Open Inventor compatibility.
#include <Inventor/SbBasic.h>
#include <Inventor/SbLinear.h>
// these list classes are defined in SbPList.h in SGI Inventor
#include <Inventor/lists/SbIntList.h>
#include <Inventor/lists/SbStringList.h>
#include <Inventor/lists/SbVec3fList.h>
// Here's the class definition of SbPList in Coin.
#include <Inventor/lists/SbPList.h>
#else // COIN_INTERNAL
#error "Do not include Inventor/SbPList.h internally (use Inventor/lists/SbPList.h)."
#endif // COIN_INTERNAL

#endif // !COIN_SBPLIST_H
