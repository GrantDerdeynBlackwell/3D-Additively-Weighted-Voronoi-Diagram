#ifndef G_POWER
#define G_POWER
#include "typedefs.h"
#include <boost/program_options/variables_map.hpp>

void power (const Model &model, Rt &T);

std::array<double, 4>
subdivide (const Rt::Vertex_handle &vh, const Rt &T,
           const boost::program_options::variables_map &vm,
           const std::set<std::string> &residues);

#endif
