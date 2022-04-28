#include "awVd_ray.h"
#include "bfs_for_vert.h"
#include "compute_bisectors.h"
#include "icosphere.h"
#include "ray_intersections.h"
#include "typedefs.h"
#include "voronoi_mesh_details.h"

#include "/home/e5-2690/builds/overlap/overlap.hpp"

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Intersections_2/Ray_2_Segment_2.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/exceptions.h>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <cstddef>
#include <deque>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>

#include <CGAL/IO/Color.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

namespace po = boost::program_options;
typedef CGAL::IO::Color Color;

void
connect (Mesh::Halfedge_index h, Icosphere &ico)
{
  while (ico.m.degree (ico.m.face (h)) > 3)
    {
      CGAL::Euler::split_face (h, ico.m.next (ico.m.next (h)), ico.m);
    }
}

NT
compute_curvature (const Hyperbola &h, const NT x, const NT y, const NT z)
{
  const NT A{ h[0] };
  const NT B{ h[1] };
  const NT C{ h[2] };
  const NT D{ h[3] };
  const NT E{ h[4] };
  const NT F{ h[5] };
  const NT G{ h[6] };
  const NT H{ h[7] };
  const NT I{ h[8] };
  const NT J{ h[9] };

  const NT dx{ G + 2. * A * x + D * y + E * z };
  const NT dy{ H + D * x + 2. * B * y + F * z };
  const NT dz{ I + E * x + F * y + 2. * C * z };

  // clang-format off
    Eigen::Matrix<NT, 4, 4> curvature_mat;
    curvature_mat << 2.*A ,    D ,    E , dx,
                        D , 2.*B ,    F , dy,
                        E ,    F , 2.*C , dz,
                       dx ,   dy ,   dz , 0.;
  // clang-format on

  // norm of the gradient to the 4th power
  NT grad_norm_4
      = std::pow ((std::pow (dx, 2) + std::pow (dy, 2) + std::pow (dz, 2)), 2);

  return curvature_mat.determinant () / grad_norm_4;
}

std::array<double, 5>
compute_volume (const Icosphere &ico, const std::set<std::string> &residues)
{
  auto vatom = ico.m.property_map<vertex_descriptor, std::set<const Atom *> > (
      "v:atom");
  assert (vatom.second);
  Mesh::Property_map<Mesh::Face_index, CGAL::IO::Color> fcolor;
  Sphere s{ vector_t{ 0., 0., 0. }, 1 };

  double vol = 0.;
  double overlap_vol = 0.;
  double surface_area = 0.;
  double interfacial_area = 0.;
  double max_K = 0.;

  for (auto face : ico.m.faces ())
    {
      bool interface_flag = false;
      for (const auto &v : ico.m.vertices_around_face (ico.m.halfedge (face)))
        {
          if (!vatom.first[v].empty ())
            {
              if (residues.find ((*vatom.first[v].begin ())->residue_name ())
                  == residues.cend ())
                {
                  if (vatom.first[v].size () == 1)
                    {
                      interface_flag = true;
                    }
                }
              for (const auto atom : vatom.first[v])
                {
                  max_K = std::max (
                      std::fabs (compute_curvature (
                          ico.h ().at ((const Atom *const)atom),
                          ico.m.point (v).x (), ico.m.point (v).y (),
                          ico.m.point (v).z ())),
                      max_K);
                }
            }
        }
      auto p0 = ico.m.point (ico.m.target (ico.m.halfedge (face)));
      auto p1
          = ico.m.point (ico.m.target (ico.m.next (ico.m.halfedge (face))));
      auto p2 = ico.m.point (
          ico.m.target (ico.m.next (ico.m.next (ico.m.halfedge (face)))));

      K::Triangle_3 tri{ p0, p1, p2 };

      surface_area += std::sqrt (tri.squared_area ());

      interfacial_area
          += interface_flag ? std::sqrt (tri.squared_area ()) : 0.;

      vector_t v0 = { 0., 0., 0. };
      vector_t v1 = { p0[0], p0[1], p0[2] };
      vector_t v2 = { p1[0], p1[1], p1[2] };
      vector_t v3 = { p2[0], p2[1], p2[2] };

      if ((v1 - v0).cross (v2 - v0).dot (v3 - v0) >= 1e-9)
        {
          Tetrahedron tet{ v0, v1, v2, v3 };
          vol += tet.volume;

          overlap_vol += tet.volume > 1e-6 ? overlap (s, tet) : 0;
        }
      else if ((v2 - v0).cross (v1 - v0).dot (v3 - v0) >= 1e-9)
        {
          Tetrahedron tet{ v0, v2, v1, v3 };
          vol += tet.volume;

          overlap_vol += tet.volume > 1e-6 ? overlap (s, tet) : 0;
        }
    }
  return { vol, overlap_vol, surface_area, interfacial_area, max_K };
}

#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
void
mesh_io_test (Icosphere ico, const std::string &fname)
{
  auto vatom = ico.m.property_map<vertex_descriptor, std::set<const Atom *> > (
      "v:atom");
  assert (vatom.second);
  Mesh::Property_map<Mesh::Face_index, CGAL::IO::Color> fcolor;

  bool fcreated;
  boost::tie (fcolor, fcreated)
      = ico.m.add_property_map<Mesh::Face_index, CGAL::IO::Color> (
          "f:color", CGAL::IO::Color ());

  for (const auto &v : ico.m.vertices ())
    {
      auto p = ico.m.point (v);
      ico.m.point (v)
          = K::Point_3 (p.x () * ico.radius () + ico.center ().x (),
                        p.y () * ico.radius () + ico.center ().y (),
                        p.z () * ico.radius () + ico.center ().z ());
    }
  for (const auto &f : ico.m.faces ())
    {
      std::set<const Atom *> intersect
          = vatom.first[ico.m.target (ico.m.halfedge (f))];
      for (const auto &v : ico.m.vertices_around_face (ico.m.halfedge (f)))
        {
          std::set<const Atom *> tmp;
          std::set_intersection (vatom.first[v].cbegin (),
                                 vatom.first[v].cend (), intersect.cbegin (),
                                 intersect.cend (),
                                 std::inserter (tmp, tmp.begin ()));

          intersect = tmp;
        }
      if (intersect.empty ())
        {
          fcolor[f] = CGAL::IO::Color{ 0, 0, 0 };
        }
      else if (intersect.size () > 1)
        {
          fcolor[f] = CGAL::IO::Color{ 0, 0, 0 };
        }
      else
        {
          auto color = intersect.empty () || intersect.size () > 1
                           ? G_cmap.end ()
                           : G_cmap.find ((*intersect.begin ())->element ());
          fcolor[f] = color != G_cmap.end () ? color->second
                                             : CGAL::IO::Color{ 0, 0, 0 };
        }
    }
  std::ofstream ofile (fname, std::ios::binary);

  CGAL::IO::write_OFF (ofile, ico.m,
                       CGAL::parameters::face_color_map (fcolor));
  ofile.close ();
  ico.m.clear ();
}
void
compute_voronoi_faces (Icosphere &icosphere)
{

  auto vcolor
      = icosphere.m.property_map<vertex_descriptor, std::set<const Atom *> > (
          "v:atom");
  assert (vcolor.second);

  // For each vertex, draw a ray and find the closest valid intersection.
  for (auto &vertex : icosphere.m.vertices ())
    {
      Ray ray{ icosphere.m.point (vertex)[0], icosphere.m.point (vertex)[1],
               icosphere.m.point (vertex)[2] };
      auto t = compute_scalar_dist (ray, icosphere.h (), icosphere.atom ());
      // Because the icosphere has already been made, just update the
      // coordinates of the vertex and the mesh will stay valid.

      Eigen::Vector3d new_coords = ray * t.first;
      icosphere.m.point (vertex)
          = K::Point_3 (new_coords[0], new_coords[1], new_coords[2]);

      if (t.second != nullptr)
        {
          vcolor.first[vertex].insert (t.second);
        }
    }
}

const std::map<int, std::string> G_omap{ { 0, "\x1B[31m" }, { 1, "\x1B[32m" },
                                         { 2, "\x1B[33m" }, { 3, "\x1B[34m" },
                                         { 4, "\x1B[35m" }, { 5, "\x1B[36m" },
                                         { 6, "\x1B[37m" } };

void
compute_voronoi_vertices (Icosphere &ico, Voronoi_map &voronoi_vertices)
{
  auto vcolor
      = ico.m.property_map<vertex_descriptor, std::set<const Atom *> > (
          "v:atom");
  assert (vcolor.second);

  auto evisited = ico.m.property_map<Mesh::Edge_index, bool> ("e:visited");
  assert (evisited.second);

  // A new property map to keep track of the visited edges when we bfs

  // A new property map to keep track of the Voronoi vertices that are
  // contained in the faces. This is better than changing the topology
  // because we can delay connecting things until we know all the faces that
  // intersect in this triangle.
  Mesh::edge_iterator eit = ico.m.edges_begin ();

  std::set<const Atom *> uni;
  while (eit != ico.m.edges_end ())
    {
      int trys = 0;
      do
        {
          trys++;
          assert (eit->is_valid ());
          uni.clear ();
          uni.insert (vcolor.first[ico.m.vertex (*eit, 0)].begin (),
                      vcolor.first[ico.m.vertex (*eit, 0)].end ());
          uni.insert (vcolor.first[ico.m.vertex (*eit, 1)].begin (),
                      vcolor.first[ico.m.vertex (*eit, 1)].end ());

          std::set<const Atom *> intersect;
          intersect.clear ();

          std::set_intersection (
              vcolor.first[ico.m.vertex (*eit, 0)].cbegin (),
              vcolor.first[ico.m.vertex (*eit, 0)].cend (),
              vcolor.first[ico.m.vertex (*eit, 1)].cbegin (),
              vcolor.first[ico.m.vertex (*eit, 1)].cend (),
              std::inserter (intersect, intersect.begin ()));

          if (uni.size () > 2 && intersect.empty ())
            {
              auto p = ico.m.point (ico.m.vertex (*eit, 0));
              auto q = ico.m.point (ico.m.vertex (*eit, 1));
              // printf ("uni.size () > 2 && intersect.empty ()\n");
              auto nh = CGAL::Euler::split_edge (ico.m.halfedge (*eit), ico.m);

              Eigen::Vector3d midpoint{ (p[0] + q[0]) / 2., (p[1] + q[1]) / 2.,
                                        (p[2] + q[2]) / 2. };

              midpoint.normalize ();

              auto t = compute_scalar_dist (midpoint, ico.h (), ico.atom ());
              ico.m.point (ico.m.target (nh))
                  = K::Point_3{ midpoint[0] * t.first, midpoint[1] * t.first,
                                midpoint[2] * t.first };

              if (t.second != nullptr)
                {
                  vcolor.first[ico.m.target (nh)].insert (t.second);
                }
              //              CGAL::draw(ico.m);
              connect (nh, ico);
              connect (ico.m.prev (ico.m.opposite (nh)), ico);
              // mesh_io_test(ico, std::to_string(ico.m.num_vertices()) +
              // ".off");
            }
          else if (intersect.empty ()
                   && !vcolor.first[ico.m.vertex (*eit, 0)].empty ()
                   && !vcolor.first[ico.m.vertex (*eit, 1)].empty ())
            {
              // printf ("intersect.empty()\n");
              bfs_for_vert (ico, *eit, vcolor.first, evisited.first,
                            voronoi_vertices);
            }
          else
            {
              break;
            }
        }
      while (uni.size () > 2);
      ++eit;
    }
}

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/polygon_mesh_processing.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

std::array<double, 5>
find_neighbors (const Model &model, const Atom &atom,
                const std::size_t subdivisions, const po::variables_map &vm,
                const std::set<std::string> &residues,
                Voronoi_map &voronoi_vertices)
{
  Hyperbola_map h;
  double r = G_atom_classifier.get_properties (atom).value ();
  for (Atom_Const_Iter it = model.atoms_begin (); it != model.atoms_end ();
       ++it)
    {
      // Filter the hyperbolas by distance. Atoms whose weighted distance is
      // more than twice the cutoff will not contribute to the cell.
      if (std::sqrt (std::pow ((it->x () - atom.x ()), 2)
                     + std::pow ((it->y () - atom.y ()), 2)
                     + std::pow ((it->z () - atom.z ()), 2))
                  - r - G_atom_classifier.get_properties (*it).value ()
              < 2. * CUTOFF
          && it->atom_serial_number () != atom.atom_serial_number ())
        {
          Hyperbola hyp = compute_bisector (atom, *it);
          h.insert ({ (const Atom *const) & *it, hyp });
        }
    }

  Icosphere ico{ atom, subdivisions, h };
  // printf ("computing faces...\n");
  compute_voronoi_faces (ico);
  // mesh (ico);
  // printf ("finished computing faces.\ncomputing vertices...\n");
  compute_voronoi_vertices (ico, voronoi_vertices);
  // mesh (ico);
  // printf ("finished computing vertices.\ncomputing edges...");
  // compute_voronoi_edges (ico);
  // printf ("finished computing edges.\n");

  if (vm.count ("mv"))
    {
      mesh_io_test (ico, "outputs/awVd_" + atom.atom_name () + "_"
                             + std::to_string (atom.atom_serial_number ())
                             + ".off");
    }
  auto v = compute_volume (ico, residues);

  h.clear ();

  ico.m.clear ();

  return { v[0] * std::pow (ico.radius (), 3),
           v[1] * std::pow (ico.radius (), 3),
           v[2] * std::pow (ico.radius (), 2),
           v[3] * std::pow (ico.radius (), 2), v[4] };
}
