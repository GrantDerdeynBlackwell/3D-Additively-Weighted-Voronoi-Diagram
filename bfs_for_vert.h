#ifndef G_BFS_FOR_VERT
#define G_BFS_FOR_VERT

#include "icosphere.h"
#include "typedefs.h"
#include <Eigen/src/Core/Matrix.h>

void bfs_for_vert (
    Icosphere &ico, Mesh::Edge_index s,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited);

#endif
