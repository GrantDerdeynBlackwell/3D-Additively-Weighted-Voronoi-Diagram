#ifndef G_COMPUTE_BISECTORS
#define G_COMPUTE_BISECTORS
#include "typedefs.h"
#include "icosphere.h"

Hyperbola const compute_bisector (const Icosphere &ico, const Atom &atom,
                                  const Atom &it);

#endif
