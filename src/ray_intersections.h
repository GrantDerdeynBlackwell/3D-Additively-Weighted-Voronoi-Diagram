#ifndef G_RAY_INTERSECTIONS
#define G_RAY_INTERSECTIONS

#include "typedefs.h"
#include "parser.h"

std::pair<double, const Atom *>
compute_scalar_dist (const Ray &ray,
                     const Hyperbola_map &h,
                     const Atom &atom);
bool
almost_equal (double x, double y, std::size_t ulp);

#endif
