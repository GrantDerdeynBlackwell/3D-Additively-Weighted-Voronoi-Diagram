#ifndef G_AWVD_RAY
#define G_AWVD_RAY

#include "icosphere.h"
#include "parser.h"
#include "typedefs.h"

#include <ESBTL/PDB.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/default.h>

void
connect (Mesh::Halfedge_index h, Icosphere &ico);

void mesh_io_test (Icosphere ico, const std::string &fname);

std::array<double, 5>
find_neighbors (const Model &model, const Atom &atom,
                const std::size_t subdivisions, const Cmd_line_options &options);

#endif
