/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_MEASUREMENTS_H
#define _LLKA_MEASUREMENTS_H

#include "llka_structure.h"

LLKA_API float LLKA_CC LLKA_measureAnglef(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c);
LLKA_API double LLKA_CC LLKA_measureAngle(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c);

LLKA_API float LLKA_CC LLKA_measureDihedralf(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c, const LLKA_Atom *d);
LLKA_API double LLKA_CC LLKA_measureDihedral(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c, const LLKA_Atom *d);

LLKA_API float LLKA_CC LLKA_measureDistancef(const LLKA_Atom *a, const LLKA_Atom *b);
LLKA_API double LLKA_CC LLKA_measureDistance(const LLKA_Atom *a, const LLKA_Atom *b);

#endif /* _LLKA_MEASUREMENTS_H */
