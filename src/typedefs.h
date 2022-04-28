#ifndef G_TYPEDEFS
#define G_TYPEDEFS

#include "ESBTL/molecular_system.h"
#include "ESBTL/xyz_utils.h"
#include <CGAL/MP_Float.h>
#include <CGAL/NT_converter.h>
#include <ESBTL/CGAL/EPIC_kernel_with_atom.h>
#include <ESBTL/PDB.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/default.h>
#include <ESBTL/occupancy_handlers.h>

#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>

#include <cstddef>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <utility>

#include <CGAL/Eigen_vector.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/STL.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Linear_algebraCd.h>

#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_svd.h>
#include <CGAL/Eigen_vector.h>

#include <stack>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>

#define CUTOFF 2.8

// awVd typedefs
// using NT = typename CGAL::MP_Float;
using NT = double;
// using K = typename CGAL::Simple_cartesian<NT>;
// using K = typename ESBTL::CGAL::EPIC_kernel_with_atom;
using K = typename CGAL::Exact_predicates_inexact_constructions_kernel;
using double_to_NT = typename CGAL::NT_converter<double, NT>;

using Items = typename ESBTL::Default_system_items;
using System = typename ESBTL::Molecular_system<Items, K::Point_3>;
using Model = typename System::Model;
using Atom = typename System::Atom;
using Atom_Iter = typename Model::Atoms_iterator;
using Atom_Const_Iter = typename Model::Atoms_const_iterator;
using T_Atom_classifier =
    typename ESBTL::Generic_classifier<ESBTL::Radius_of_atom<NT, Atom> >;

using Hyperbola = typename std::array<NT, 10>;
using Hyperbola_map = typename std::map<const Atom *const, const Hyperbola>;

using Ray = typename Eigen::Vector3d;

using Index = typename std::size_t;

using Accept_none_occupancy_policy =
    typename ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >;

const T_Atom_classifier G_atom_classifier ("./data/bondi_classifier.txt");
const double_to_NT G_double_to_NT;

using Voronoi_map = typename std::map<const std::array<const Atom *, 4>,
                                      const std::array<const double, 4> >;

// using K = typename CGAL::Simple_cartesian<double>;
using Mesh = typename CGAL::Surface_mesh<K::Point_3>;
using vertex_descriptor = typename Mesh::Vertex_index;
using face_descriptor = typename Mesh::Face_index;

using EK = typename CGAL::Exact_predicates_exact_constructions_kernel;
using Weighted_point = typename EK::Weighted_point_3;
using Vb0 = typename CGAL::Regular_triangulation_vertex_base_3<EK>;
using Cb0 = typename CGAL::Regular_triangulation_cell_base_3<EK>;

using Vb =
    typename CGAL::Triangulation_vertex_base_with_info_3<Atom *, EK, Vb0>;
using Vb_vvert = typename CGAL::Triangulation_vertex_base_with_info_3<
    std::array<const Atom *, 4>, EK, Vb0>;
using Tds = typename CGAL::Triangulation_data_structure_3<Vb, Cb0>;
using Tds_vvert = typename CGAL::Triangulation_data_structure_3<Vb_vvert, Cb0>;
using Rt = typename CGAL::Regular_triangulation_3<EK, Tds>;
using SubTri = typename CGAL::Regular_triangulation_3<EK, Tds_vvert>;

typedef CGAL::IO::Color Color;
const std::map<const std::string, Color> G_cmap
    = { { "H", Color (255, 255, 255) },  { "He", Color (217, 255, 255) },
        { "Li", Color (204, 128, 255) }, { "Be", Color (194, 255, 0) },
        { "B", Color (255, 181, 181) },  { "C", Color (144, 144, 144) },
        { "N", Color (48, 80, 248) },    { "O", Color (255, 13, 13) },
        { "F", Color (144, 224, 80) },   { "Ne", Color (179, 227, 245) },
        { "Na", Color (171, 92, 242) },  { "Mg", Color (138, 255, 0) },
        { "Al", Color (191, 166, 166) }, { "Si", Color (240, 200, 160) },
        { "P", Color (255, 128, 0) },    { "S", Color (255, 255, 48) },
        { "Cl", Color (31, 240, 31) },   { "Ar", Color (128, 209, 227) },
        { "K", Color (143, 64, 212) },   { "Ca", Color (61, 255, 0) },
        { "Sc", Color (230, 230, 230) }, { "Ti", Color (191, 194, 199) },
        { "V", Color (166, 166, 171) },  { "Cr", Color (138, 153, 199) },
        { "Mn", Color (156, 122, 199) }, { "Fe", Color (224, 102, 51) },
        { "Co", Color (240, 144, 160) }, { "Ni", Color (80, 208, 80) },
        { "Cu", Color (200, 128, 51) },  { "Zn", Color (125, 128, 176) },
        { "Ga", Color (194, 143, 143) }, { "Ge", Color (102, 143, 143) },
        { "As", Color (189, 128, 227) }, { "Se", Color (255, 161, 0) },
        { "Br", Color (166, 41, 41) },   { "Kr", Color (92, 184, 209) },
        { "Rb", Color (112, 46, 176) },  { "Sr", Color (0, 255, 0) },
        { "Y", Color (148, 255, 255) },  { "Zr", Color (148, 224, 224) },
        { "Nb", Color (115, 194, 201) }, { "Mo", Color (84, 181, 181) },
        { "Tc", Color (59, 158, 158) },  { "Ru", Color (36, 143, 143) },
        { "Rh", Color (10, 125, 140) },  { "Pd", Color (0, 105, 133) },
        { "Ag", Color (192, 192, 192) }, { "Cd", Color (255, 217, 143) },
        { "In", Color (166, 117, 115) }, { "Sn", Color (102, 128, 128) },
        { "Sb", Color (158, 99, 181) },  { "Te", Color (212, 122, 0) },
        { "I", Color (148, 0, 148) },    { "Xe", Color (66, 158, 176) },
        { "Cs", Color (87, 23, 143) },   { "Ba", Color (0, 201, 0) },
        { "La", Color (112, 212, 255) }, { "Ce", Color (255, 255, 199) },
        { "Pr", Color (217, 255, 199) }, { "Nd", Color (199, 255, 199) },
        { "Pm", Color (163, 255, 199) }, { "Sm", Color (143, 255, 199) },
        { "Eu", Color (97, 255, 199) },  { "Gd", Color (69, 255, 199) },
        { "Tb", Color (48, 255, 199) },  { "Dy", Color (31, 255, 199) },
        { "Ho", Color (0, 255, 156) },   { "Er", Color (0, 230, 117) },
        { "Tm", Color (0, 212, 82) },    { "Yb", Color (0, 191, 56) },
        { "Lu", Color (0, 171, 36) },    { "Hf", Color (77, 194, 255) },
        { "Ta", Color (77, 166, 255) },  { "W", Color (33, 148, 214) },
        { "Re", Color (38, 125, 171) },  { "Os", Color (38, 102, 150) },
        { "Ir", Color (23, 84, 135) },   { "Pt", Color (208, 208, 224) },
        { "Au", Color (255, 209, 35) },  { "Hg", Color (184, 184, 208) },
        { "Tl", Color (166, 84, 77) },   { "Pb", Color (87, 89, 97) },
        { "Bi", Color (158, 79, 181) },  { "Po", Color (171, 92, 0) },
        { "At", Color (117, 79, 69) },   { "Rn", Color (66, 130, 150) },
        { "Fr", Color (66, 0, 102) },    { "Ra", Color (0, 125, 0) },
        { "Ac", Color (112, 171, 250) }, { "Th", Color (0, 186, 255) },
        { "Pa", Color (0, 161, 255) },   { "U", Color (0, 143, 255) },
        { "Np", Color (0, 128, 255) },   { "Pu", Color (0, 107, 255) },
        { "Am", Color (84, 92, 242) },   { "Cm", Color (120, 92, 227) },
        { "Bk", Color (138, 79, 227) },  { "Cf", Color (161, 54, 212) },
        { "Es", Color (179, 31, 212) },  { "Fm", Color (179, 31, 186) },
        { "Md", Color (179, 13, 166) },  { "No", Color (189, 13, 135) },
        { "Lr", Color (199, 0, 102) },   { "Rf", Color (204, 0, 89) },
        { "Db", Color (209, 0, 79) },    { "Sg", Color (217, 0, 69) },
        { "Bh", Color (224, 0, 56) },    { "Hs", Color (230, 0, 46) },
        { "Mt", Color (235, 0, 38) } };
#endif
