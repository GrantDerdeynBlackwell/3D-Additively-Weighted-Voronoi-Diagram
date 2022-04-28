#ifndef G_VORONOI_MESH_DETAILS
#define G_VORONOI_MESH_DETAILS

#include "icosphere.h"
#include "typedefs.h"
#include <Eigen/src/Core/Matrix.h>

std::vector<Eigen::Matrix<NT, 3, 1> >
compute_voronoi_vertex_rigorous (
    const Icosphere &ico, const std::array<const Atom *, 3> atoms);

std::list<K::Point_3>
find_voronoi_curves_along_mesh_edge (
    const std::array<const Atom *, 2> &curve_colors,
    Icosphere &ico, const Mesh::Edge_index &edge);

#endif
