#ifndef G_AWVD_RAY
#define G_AWVD_RAY

#include <ESBTL/PDB.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/default.h>

#include "icosphere.h"
#include "typedefs.h"

void connect (Mesh::Halfedge_index h, Icosphere &ico);

void mesh_io_test (Icosphere ico, const std::string &fname);

std::array<double, 5>
find_neighbors (const Model &model, const Atom &atom,
                const std::size_t subdivisions,
                const boost::program_options::variables_map &vm,
                const std::set<std::string> &residues,
                Voronoi_map &voronoi_vertices,
                const std::map<const Atom *, double> &radii_map);

#endif
