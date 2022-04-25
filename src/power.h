#ifndef G_POWER
#define G_POWER
#include "parser.h"
#include "typedefs.h"

void power (const Model &model, Rt &T);

std::array<double, 4> subdivide (const Rt::Vertex_handle &vh, const Rt &T,
                                 Cmd_line_options &options);

#endif
