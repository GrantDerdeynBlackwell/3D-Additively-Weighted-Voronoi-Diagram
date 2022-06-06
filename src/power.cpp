#include "power.h"

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/PLY.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Lazy.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_3_to_face_graph.h>
#include <CGAL/exceptions.h>
#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/spatial_sort.h>
#include <ESBTL/PDB.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/default.h>

#include <algorithm>
#include <boost/graph/properties.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cassert>
#include <memory>
#include <string>
#include <vector>

#include "overlap.hpp"
#include "ESBTL/occupancy_handlers.h"
#include "typedefs.h"

namespace po = boost::program_options;

struct WeightedLessX
{
  bool
  operator() (const Weighted_point &p, const Weighted_point &q) const
  {
    return p.x () < q.x ();
  }
};
struct WeightedLessY
{
  bool
  operator() (const Weighted_point &p, const Weighted_point &q) const
  {
    return p.y () < q.y ();
  }
};
struct WeightedLessZ
{
  bool
  operator() (const Weighted_point &p, const Weighted_point &q) const
  {
    return p.z () < q.z ();
  }
};

struct Weighted_SpatialSortingTraits
{
  typedef Weighted_point Point_3;
  typedef WeightedLessX Less_x_3;
  typedef WeightedLessY Less_y_3;
  typedef WeightedLessY Less_z_3;
  Less_x_3
  less_x_3_object () const
  {
    return Less_x_3 ();
  }
  Less_y_3
  less_y_3_object () const
  {
    return Less_y_3 ();
  }
  Less_z_3
  less_z_3_object () const
  {
    return Less_z_3 ();
  }
};

template <class Triangulation_3, class FG>
typename boost::graph_traits<FG>::vertex_descriptor
link_to_fg_with_color (const Triangulation_3 &t,
                       typename Triangulation_3::Vertex_handle vh, FG &fg,
                       const Atom *atom, const std::string &fname,
                       bool no_infinite_faces = true)
{
  typedef typename Triangulation_3::Cell_handle Cell_handle;
  typedef typename Triangulation_3::Vertex_handle Vertex_handle;
  typedef
      typename boost::graph_traits<FG>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FG>::face_descriptor face_descriptor;
  typedef CGAL::dynamic_face_property_t<CGAL::IO::Color> fcolor_tag;
  typedef typename boost::property_map<FG, fcolor_tag>::type face_id_map;

  clear (fg);

  fg.clear ();

  vertex_descriptor inf;
  vertex_descriptor nullvertex = boost::graph_traits<FG>::null_vertex ();

  typedef boost::unordered_map<Vertex_handle, vertex_descriptor> Vertex_map;
  typedef boost::unordered_map<vertex_descriptor, Vertex_handle> RVertex_map;
  Vertex_map vertex_map;
  RVertex_map rvertex_map;

  std::vector<Cell_handle> cells;
  t.incident_cells (t.infinite_vertex (), std::back_inserter (cells));
  std::array<vertex_descriptor, 3> face;

  typename boost::property_map<FG, CGAL::vertex_point_t>::type vpm
      = get (CGAL::vertex_point, fg);

  face_id_map fcolor = get (fcolor_tag (), fg);

  for (Cell_handle ch : cells)
    {
      bool infinite_face = false;
      int vhi = ch->index (vh);
      for (int i = 0; i < 3; i++)
        {
          int j = Triangulation_3::vertex_triple_index (vhi, i);
          Vertex_handle vhj = ch->vertex (j);
          if (no_infinite_faces && t.is_infinite (vhj))
            {
              infinite_face = true;
            }
          else
            {
              std::pair<typename Vertex_map::iterator, bool> res
                  = vertex_map.insert (std::make_pair (vhj, nullvertex));
              if (res.second)
                {
                  res.first->second = add_vertex (fg);
                  put (vpm, res.first->second, vhj->point ().point ());
                  if (t.is_infinite (vhj))
                    {
                      inf = res.first->second;
                    }
                  rvertex_map.insert (std::make_pair (res.first->second, vhj));
                }
              face[i] = res.first->second;
            }
        }
      if (!infinite_face)
        {
          std::set<const Atom *> intersect{
            rvertex_map[face[0]]->info ().cbegin (),
            rvertex_map[face[0]]->info ().cend ()
          };
          for (const auto &v : face)
            {
              std::set<const Atom *> tmp;
              std::set_intersection (rvertex_map[v]->info ().cbegin (),
                                     rvertex_map[v]->info ().cend (),
                                     intersect.cbegin (), intersect.cend (),
                                     std::inserter (tmp, tmp.begin ()));

              intersect = tmp;
            }

          intersect.erase (atom);

	  if(intersect.empty ())
	  {
	    throw std::runtime_error("This input has one or more unbounded power cells. It is not currently possible to mesh unbounded power cells.");
	  }
          auto color = G_cmap.find ((*intersect.cbegin ())->element ());
          auto new_face = CGAL::Euler::add_face (face, fg);
          put (fcolor, new_face,
               color != G_cmap.cend () ? color->second
                                       : CGAL::IO::Color{ 0, 0, 0 });
        }
    }
  std::ofstream ofile (fname, std::ios::binary);
  CGAL::IO::write_OFF (ofile, fg, CGAL::parameters::face_color_map (fcolor));
  return inf;
}

void
mesh (const SubTri &T, const std::string &fname, const Atom *atom)
{
  CGAL::Surface_mesh<EK::Point_3> m;
  link_to_fg_with_color (T, T.infinite_vertex (), m, atom, fname);
}

void
power (const Model &model, Rt &T,
       const std::map<const Atom *, double> &radii_map)
{
  std::vector<std::pair<Weighted_point, Atom *> > points;
  for (Atom_Const_Iter it = model.atoms_begin (); it != model.atoms_end ();
       ++it)
    {
      points.push_back (std::make_pair (
          Weighted_point (EK::Point_3{ it->x (), it->y (), it->z () },
                          std::pow (radii_map.at (&*it), 2)),
          (Atom *)&*it));
    }
  T.insert (points.begin (), points.end ());

  assert (T.is_valid ());
  assert (T.dimension () == 3);
}

std::array<double, 4>
subdivide (const Rt::Vertex_handle &vh, const Rt &T,
           const po::variables_map &vm, const std::set<std::string> &residues)
{
  double vol = 0.;
  double overlap_vol = 0.;
  double area = 0.;
  double interfacial_area = 0.;
  std::vector<Rt::Cell_handle> cells;
  T.incident_cells (
      vh, std::back_insert_iterator<std::vector<Rt::Cell_handle> > (cells));

  std::vector<std::pair<Weighted_point, std::array<const Atom *, 4> > >
      poly_points;
  for (auto cell : cells)
    {
      std::array<const Atom *, 4> colors{ cell->vertex (0)->info (),
                                          cell->vertex (1)->info (),
                                          cell->vertex (2)->info (),
                                          cell->vertex (3)->info () };

      std::sort (colors.begin (), colors.end ());
      if (!T.is_infinite (cell))
        {
          poly_points.push_back ({ { T.dual (cell), 1. }, colors });
        }
    }
  SubTri subtri{ poly_points.begin (), poly_points.end () };

  // Don't forget to sqrt the weight to get the radius
  Sphere s{ vector_t{ CGAL::to_double (vh->point ().x ()),
                      CGAL::to_double (vh->point ().y ()),
                      CGAL::to_double (vh->point ().z ()) },
            std::sqrt (CGAL::to_double (vh->point ().weight ())) };

  for (const auto &facet : subtri.finite_facets ())
    {
      std::vector<int> points{ 0, 1, 2, 3 };
      points.erase (std::remove (points.begin (), points.end (), facet.second),
                    points.end ());
      const auto &v0 = facet.first->vertex (facet.second);
      const auto &v1 = facet.first->vertex (points[0]);
      const auto &v2 = facet.first->vertex (points[1]);
      const auto &v3 = facet.first->vertex (points[2]);

      // The facet lies on the surface of the Voronoi cell if it is incident to
      // the infinite cell.

      CGAL::Object test = subtri.dual (facet);
      if (const EK::Ray_3 *r = CGAL::object_cast<EK::Ray_3> (&test))
        {
          CGAL::Triangle_3<EK> tri (v1->point ().point (),
                                    v2->point ().point (),
                                    v3->point ().point ());

          std::set<const Atom *> cell_color;
          std::set_intersection (
              v1->info ().cbegin (), v1->info ().cend (),
              v2->info ().cbegin (), v2->info ().cend (),
              std::inserter (cell_color, cell_color.begin ()));
          std::set_intersection (
              v1->info ().cbegin (), v1->info ().cend (),
              v3->info ().cbegin (), v3->info ().cend (),
              std::inserter (cell_color, cell_color.begin ()));
          bool interface_flag = false;
          for (const auto &elem : cell_color)
            {
              if (residues.find (elem->residue_name ()) == residues.cend ())
                {
                  interface_flag = true;
                }
            }
          area += std::sqrt (CGAL::to_double (tri.squared_area ()));
          interfacial_area
              += interface_flag
                     ? std::sqrt (CGAL::to_double (tri.squared_area ()))
                     : 0.;

          vector_t tet0 = { CGAL::to_double (vh->point ().x ()),
                            CGAL::to_double (vh->point ().y ()),
                            CGAL::to_double (vh->point ().z ()) };
          vector_t tet1 = { CGAL::to_double (v1->point ().x ()),
                            CGAL::to_double (v1->point ().y ()),
                            CGAL::to_double (v1->point ().z ()) };
          vector_t tet2 = { CGAL::to_double (v2->point ().x ()),
                            CGAL::to_double (v2->point ().y ()),
                            CGAL::to_double (v2->point ().z ()) };
          vector_t tet3 = { CGAL::to_double (v3->point ().x ()),
                            CGAL::to_double (v3->point ().y ()),
                            CGAL::to_double (v3->point ().z ()) };

          if ((tet1 - tet0).cross (tet2 - tet0).dot (tet3 - tet0) >= 1e-9)
            {
              Tetrahedron tet{ tet0, tet1, tet2, tet3 };
              vol += tet.volume;

              overlap_vol += tet.volume > 1e-9 ? overlap (s, tet) : 0;
            }
          else if ((tet2 - tet0).cross (tet1 - tet0).dot (tet3 - tet0) >= 1e-9)
            {
              Tetrahedron tet{ tet0, tet2, tet1, tet3 };
              vol += tet.volume;

              overlap_vol += tet.volume > 1e-9 ? overlap (s, tet) : 0;
            }
        }
    }
  if (vm.count ("mp"))
    {
      std::string fname{ "outputs/power_" + vh->info ()->atom_name () + "_"
                         + std::to_string (vh->info ()->atom_serial_number ())
                         + ".off" };
      mesh (subtri, fname, vh->info ());
    }
  return { vol, overlap_vol, area, interfacial_area };
}
