/* vim: set sw=4 ts=4 sts=4 expandtab : */



#include <llka_measurements.h>

#include <util/geometry.h>



float LLKA_CC LLKA_measureAnglef(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c)
{
    return LLKAInternal::angle<float>(a->coords, b->coords, c->coords);
}


double LLKA_CC LLKA_measureAngle(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c)
{
    return LLKAInternal::angle<double>(a->coords, b->coords, c->coords);
}


float LLKA_CC LLKA_measureDihedralf(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c, const LLKA_Atom *d)
{
    return LLKAInternal::dihedralAngle<float>(a->coords, b->coords, c->coords, d->coords);
}


double LLKA_CC LLKA_measureDihedral(const LLKA_Atom *a, const LLKA_Atom *b, const LLKA_Atom *c, const LLKA_Atom *d)
{
    return LLKAInternal::dihedralAngle<double>(a->coords, b->coords, c->coords, d->coords);
}


float LLKA_CC LLKA_measureDistancef(const LLKA_Atom *a, const LLKA_Atom *b)
{
    return LLKAInternal::spatialDistance<float>(*a, *b);
}


double LLKA_CC LLKA_measureDistance(const LLKA_Atom *a, const LLKA_Atom *b)
{
    return LLKAInternal::spatialDistance<double>(*a, *b);
}
