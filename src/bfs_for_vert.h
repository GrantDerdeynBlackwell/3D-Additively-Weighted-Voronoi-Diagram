#ifndef G_BFS_FOR_VERT
#define G_BFS_FOR_VERT

#include "icosphere.h"
#include "typedefs.h"
#include <Eigen/src/Core/Matrix.h>

bool
handle_voronoi_vertex (
    Icosphere &ico, const std::array<const Atom *, 2> &curve_colors,
    const Atom *new_color,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    Mesh::Halfedge_index h, int &trys, Voronoi_map &voronoi_vertices);

void bfs_for_vert (
    Icosphere &ico, Mesh::Edge_index s,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    Voronoi_map &voronoi_vertices);

#endif
