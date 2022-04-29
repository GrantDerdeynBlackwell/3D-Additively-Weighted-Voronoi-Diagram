#ifndef G_RAY_INTERSECTIONS
#define G_RAY_INTERSECTIONS

#include "icosphere.h"
#include "typedefs.h"

std::pair<double, const Atom *>
compute_scalar_dist (const Ray &ray,
                     const Icosphere &ico);
bool
almost_equal (double x, double y, std::size_t ulp);

#endif
