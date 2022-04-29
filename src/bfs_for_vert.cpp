#include "bfs_for_vert.h"
#include "awVd_ray.h"
#include "icosphere.h"
#include "ray_intersections.h"
#include "typedefs.h"
#include "voronoi_mesh_details.h"

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Intersections_3/Ray_3_Segment_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/centroid.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/squared_distance_2_1.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/VectorBlock.h>
#include <algorithm>
#include <boost/optional/optional_io.hpp>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <queue>
#include <stdexcept>
#include <string>

using namespace std;

void
store_voronoi_vertex (const Icosphere &ico, Voronoi_map &voronoi_vertices,
                      const std::array<const Atom *, 4> &v,
                      const Eigen::Vector3d &v_vert)
{
  auto x = v_vert[0] * ico.radius () + ico.center ()[0];
  auto y = v_vert[1] * ico.radius () + ico.center ()[1];
  auto z = v_vert[2] * ico.radius () + ico.center ()[2];
  auto ok = voronoi_vertices.insert ({ v, { { x, y, z, 0. } } });

  if (!ok.second)
    {
      // In some situations, the same atoms can have 2 Voronoi
      // vertices. I use a very lax interval (1e-5) because the
      // vertices are the results of linear algebra calculations
      // and they are multiplied.
      auto &points = ok.first->second;
      bool is_present = false;
      for (const auto &point : points)
        {
          if (fabs (point[0] - x) < 1e-5 && fabs (point[1] - y) < 1e-5
              && fabs (point[2] - z) < 1e-5)
            {
              is_present = true;
            }
        }
      if (!is_present)
        {
          points.push_back ({ x, y, z, 0. });
        }
    }
}

void
connect_and_visit (Mesh::Halfedge_index h, Icosphere &ico,
                   Mesh::Property_map<Mesh::Edge_index, bool> &evisited)
{
  while (ico.m.degree (ico.m.face (h)) > 3)
    {
      auto nh
          = CGAL::Euler::split_face (h, ico.m.next (ico.m.next (h)), ico.m);
      evisited[ico.m.edge (nh)] = true;
    }
}
Mesh::Halfedge_index
face_has_v_vert (
    Icosphere &ico,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    const std::array<const Atom *, 2> &curve_colors, Mesh::Halfedge_index h)
{
  CGAL::Halfedge_around_face_circulator<Mesh> hbegin (h, ico.m),
      hdone (hbegin);
  ++(hbegin);
  do
    {
      std::set<const Atom *> intersect;

      std::set_intersection (curve_colors.cbegin (), curve_colors.cend (),
                             vcolor[ico.m.target (*hbegin)].cbegin (),
                             vcolor[ico.m.target (*hbegin)].cend (),
                             std::inserter (intersect, intersect.begin ()));

      if (intersect.size () >= 2)
        {
          return *hbegin;
        }
      intersect.clear ();
      hbegin++;
    }
  while (hbegin != hdone);

  return Mesh::null_halfedge ();
}

std::pair<NT, const Atom *>
get_t_if_exists (
    Icosphere &ico, const Mesh::Halfedge_index h,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    const Eigen::Matrix<NT, 3, 1> &v_normal, NT v_norm,
    const std::array<const Atom *, 3> atoms, Voronoi_map &voronoi_vertices)
{
  K::Plane_3 plane{ -v_normal[0], -v_normal[1], -v_normal[2], v_norm };
  CGAL::Polygon_2<K> poly;
  NT v_x = v_normal.x () * v_norm;
  NT v_y = v_normal.y () * v_norm;
  NT v_z = v_normal.z () * v_norm;
  for (auto h_edge : ico.m.halfedges_around_face (h))
    {
      auto p = ico.m.point (ico.m.target (h_edge));
      if (*vcolor[ico.m.target (h_edge)].cbegin () == nullptr)
        {
          return { 0., nullptr };
        }

      NT p_x = p.x ();
      NT p_y = p.y ();
      NT p_z = p.z ();

      bool eq_x = fabs (p_x - v_x) < 1e-6;
      bool eq_y = fabs (p_y - v_y) < 1e-6;
      bool eq_z = fabs (p_z - v_z) < 1e-6;

      if (eq_x && eq_y && eq_z)
        {
          vcolor[ico.m.target (h_edge)].clear ();
          vcolor[ico.m.target (h_edge)].insert (atoms.cbegin (),
                                                atoms.cend ());

          std::array<const Atom *, 4> v{ atoms[0], atoms[1], atoms[2],
                                         &ico.atom () };
          std::sort (v.begin (), v.end ());
          store_voronoi_vertex (ico, voronoi_vertices, v, { v_x, v_y, v_z });
          for (auto hat : ico.m.halfedges_around_target (h_edge))
            {
              connect_and_visit (hat, ico, evisited);
            }

          return { 0., nullptr };
        }

      poly.push_back (plane.to_2d (p));
      //      evisited[ico.m.edge (h_edge)] = true;
    }
  if (!poly.is_simple ())
    {
      return { 0., nullptr };
    }
  auto ray_plane = CGAL::intersection (
      K::Ray_3{ K::Point_3{ 0., 0., 0. },
                K::Point_3{ v_normal[0], v_normal[1], v_normal[2] } },
      plane);

  if (!ray_plane.is_initialized ())
    {
      return { 0., nullptr };
    }

  const K::Point_3 *p = boost::get<K::Point_3> (&*ray_plane);
  if (p)
    {
      auto two_d_p = plane.to_2d (*p);
      bool intersection = poly.has_on_bounded_side (two_d_p);
      if (!intersection)
        {
          return { 0., nullptr };
        }
      return compute_scalar_dist (Ray{ v_normal[0], v_normal[1], v_normal[2] },
                                  ico);
    }
  return { 0., nullptr };
}

Mesh::Halfedge_index
try_to_add_vert (
    const Eigen::Matrix<NT, 3, 1> &v_vert, Icosphere &ico,
    const std::array<const Atom *, 3> atoms,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    Mesh::Halfedge_index h, std::vector<const Atom *> &backups,
    Voronoi_map &voronoi_vertices)
{
  NT v_norm = v_vert.norm ();
  Eigen::Matrix<NT, 3, 1> v_normal = v_vert.normalized ();

  auto t = get_t_if_exists (ico, h, vcolor, evisited, v_normal, v_norm, atoms,
                            voronoi_vertices);

  if (t.second == nullptr)
    {
      return Mesh::null_halfedge ();
    }
  else if (std::find (atoms.cbegin (), atoms.cend (), t.second)
           == atoms.cend ())
    {
      backups.push_back (t.second);
      return Mesh::null_halfedge ();
    }
  else
    {
      // It is possible that there is another valid Voronoi vertex in
      // this triangle.
      auto nv = CGAL::Euler::add_center_vertex (h, ico.m);
      ico.m.point (ico.m.target (nv))
          = K::Point_3{ v_normal[0] * t.first, v_normal[1] * t.first,
                        v_normal[2] * t.first };
      vcolor[ico.m.target (nv)].insert (atoms.cbegin (), atoms.cend ());
      return nv;
    }
}

bool
handle_voronoi_vertex (
    Icosphere &ico, const std::array<const Atom *, 2> &curve_colors,
    const Atom *new_color,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    Mesh::Halfedge_index h, int &trys, Voronoi_map &voronoi_vertices)
{

  std::array<const Atom *, 3> atoms{ curve_colors[0], curve_colors[1],
                                     new_color };
  std::sort (atoms.begin (), atoms.end ());
  bool ok = false;
  // Something is wrong. Just return and deal with it later.
  if (trys > 1000)
    {
      return ok;
    }

  auto v_verts = compute_voronoi_vertex_rigorous (ico, atoms);

  std::vector<const Atom *> backups;

  // Create a range of face indices to be evaluated, starting with the current
  // face.
  std::set<Mesh::Face_index> f_range{ ico.m.face (h) };

  while (!v_verts.empty ())
    {
      auto v_vert = v_verts.back ();
      v_verts.pop_back ();
      for (const auto feval : f_range)
        {
          const auto heval = ico.m.halfedge (feval);
          auto nv = try_to_add_vert (v_vert, ico, atoms, vcolor, evisited,
                                     heval, backups, voronoi_vertices);
          if (nv != Mesh::null_halfedge ())
            {
              ok = true;
              std::array<const Atom *, 4> v{ curve_colors[0], curve_colors[1],
                                             new_color, &ico.atom () };
              std::sort (v.begin (), v.end ());

              store_voronoi_vertex (ico, voronoi_vertices, v, v_vert);

              auto faces_around_nv = ico.m.faces_around_target (nv);
              // If a Voronoi vertex was added, there are now new
              // triangles. We add them to f_range so they can be evaluated
              // as well.
              f_range.insert (faces_around_nv.begin (),
                              faces_around_nv.end ());
              break;
            }
        }
    }

  for (const auto &backup : backups)
    {
      for (const auto feval : f_range)
        {
          handle_voronoi_vertex (ico, curve_colors, backup, vcolor, evisited,
                                 ico.m.halfedge (feval), ++trys,
                                 voronoi_vertices);
        }
    }

  return ok;
}

void
add_voronoi_curve (
    Icosphere &ico, std::list<Mesh::Halfedge_index> &queue,
    Mesh::Halfedge_index h_int,
    const std::array<const Atom *, 2> &curve_colors,
    const K::Point_3 &v_curve_int,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited)
{
  auto nh = CGAL::Euler::split_edge (h_int, ico.m);
  // don't mark this as visited because we will start here for our next face

  ico.m.point (ico.m.target (nh))
      = K::Point_3{ v_curve_int[0], v_curve_int[1], v_curve_int[2] };
  vcolor[ico.m.target (nh)].insert (curve_colors.cbegin (),
                                    curve_colors.cend ());

  queue.push_back (ico.m.opposite (h_int));

  connect_and_visit (nh, ico, evisited);
}

void
trace_out_edge (
    Icosphere &ico, const std::array<const Atom *, 2> &curve_colors,
    Mesh::Halfedge_index s,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    Voronoi_map &voronoi_vertices)
{
  for (const auto edge : ico.m.edges ())
    {
      evisited[edge] = false;
    }
  list<Mesh::Halfedge_index> queue;
  queue.push_back (s);
  evisited[ico.m.edge (s)] = false;
  evisited[ico.m.edge (ico.m.next (s))] = false;
  while (!queue.empty ())
    {
      const auto h = queue.front ();
      queue.pop_front ();

      if (evisited[ico.m.edge (h)])
        {
          continue;
        }

      evisited[ico.m.edge (h)] = true;
      evisited[ico.m.edge (ico.m.next (h))] = true;

      auto test_h = face_has_v_vert (ico, vcolor, curve_colors, h);

      if (test_h != Mesh::null_halfedge ())
        {
          connect (test_h, ico);
          continue;
        }

      const std::vector<Mesh::Halfedge_index> hs_to_eval{
        ico.m.halfedges_around_face (h).begin (),
        ico.m.halfedges_around_face (h).end ()
      };

      std::set<const Atom *> diff;

      for (auto &heval : hs_to_eval)
        {
          assert (h.is_valid ());
          if (evisited[ico.m.edge (heval)] == true
              || vcolor[ico.m.source (heval)].find (nullptr)
                     != vcolor[ico.m.source (heval)].end ()
              || vcolor[ico.m.target (heval)].find (nullptr)
                     != vcolor[ico.m.target (heval)].end ())
            {
              continue;
            }

          auto v_curve_ints = find_voronoi_curves_along_mesh_edge (
              curve_colors, ico, ico.m.edge (heval));

          std::pair<double, const Atom *> t{ 0., nullptr };

          Mesh::Halfedge_index nh = Mesh::null_halfedge ();

          for (const auto &v_curve_int : v_curve_ints)
            {
              if (v_curve_int == K::Point_3{ 0., 0., 0. })
                {
                  continue;
                }

              auto vci_norm
                  = Ray{ v_curve_int[0], v_curve_int[1], v_curve_int[2] }
                        .normalized ();
              t = compute_scalar_dist (vci_norm, ico);

              // If this curve is covered by another face, then that means
              // the vertex must lie in this triangle.
              if (std::find (curve_colors.cbegin (), curve_colors.cend (),
                             t.second)
                  == curve_colors.cend ())
                {
                  int trys = 0;

                  if (t.second == nullptr)
                    {
                      connect_and_visit (h, ico, evisited);
                      continue;
                    }
                  handle_voronoi_vertex (ico, curve_colors, t.second, vcolor,
                                         evisited, h, trys, voronoi_vertices);
                }
              // This means the curve leaves the face. It may form a loop or
              // hit the Voronoi vertex later.
              else
                {
                  auto h_int = heval;
                  if (nh != Mesh::null_halfedge ())
                    {
                      // This edge was already intersected once. The next
                      // intersection could either lie on heval or
                      // ico.m.prev (ico.m.opposite (h_int).
                      K::Plane_3 plane{ ico.m.point (ico.m.source (nh)),
                                        ico.m.point (ico.m.target (nh)),
                                        K::Point_3{ 0., 0., 0. } };

                      K::Ray_2 ray{ plane.to_2d (K::Point_3{ 0., 0., 0. }),
                                    plane.to_2d (v_curve_int) };
                      K::Segment_2 seg{
                        plane.to_2d (ico.m.point (ico.m.source (nh))),
                        plane.to_2d (ico.m.point (ico.m.target (nh)))
                      };
                      auto check_int = CGAL::intersection (ray, seg);
                      h_int = check_int ? nh : heval;
                    }
                  add_voronoi_curve (ico, queue, h_int, curve_colors,
                                     v_curve_int, vcolor, evisited);
                  if (nh != Mesh::null_halfedge ())
                    {
                      // Because this edge was already encountered, the last
                      // member in the queue is either incorrect or is a
                      // duplicate
                      connect_and_visit (ico.m.prev (ico.m.opposite (h_int)),
                                         ico, evisited);
                      queue.pop_back ();
                    }
                  nh = ico.m.prev (ico.m.opposite (h_int));
                }
            }
        }
      connect_and_visit (h, ico, evisited);
    }
}

void
bfs_for_vert (
    Icosphere &ico, Mesh::Edge_index s,
    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > &vcolor,
    Mesh::Property_map<Mesh::Edge_index, bool> &evisited,
    Voronoi_map &voronoi_vertices)
{
  for (const auto edge : ico.m.edges ())
    {
      evisited[edge] = false;
    }

  std::set<const Atom *> colors_set;
  colors_set.insert (vcolor[ico.m.vertex (s, 0)].begin (),
                     vcolor[ico.m.vertex (s, 0)].end ());
  colors_set.insert (vcolor[ico.m.vertex (s, 1)].begin (),
                     vcolor[ico.m.vertex (s, 1)].end ());

  assert (colors_set.size () == 2);
  auto h = ico.m.halfedge (s);

  const auto v1_original = vcolor[ico.m.vertex (s, 0)];
  const auto v2_original = vcolor[ico.m.vertex (s, 1)];

  for (const auto &color1 : v1_original)
    {
      for (const auto &color2 : v2_original)
        {
          std::array<const Atom *, 2> curve_colors{ color1, color2 };
          std::sort (curve_colors.begin (), curve_colors.end ());

          if (curve_colors[0] == nullptr || curve_colors[1] == nullptr)
            {
              continue;
            }

          auto v_curve_ints
              = find_voronoi_curves_along_mesh_edge (curve_colors, ico, s);
          // Add the voronoi curve to the first edge
          for (const auto &v_curve_int : v_curve_ints)
            {
              auto int_norm
                  = Ray{ v_curve_int[0], v_curve_int[1], v_curve_int[2] }
                        .normalized ();
              auto t = compute_scalar_dist (int_norm, ico);

              auto nh = CGAL::Euler::split_edge (ico.m.halfedge (s, 0), ico.m);
              // Check that the edge does not miss another face
              if (std::find (curve_colors.begin (), curve_colors.end (),
                             t.second)
                  == curve_colors.end ())
                {
                  ico.m.point (ico.m.target (nh))
                      = K::Point_3{ int_norm[0] * t.first,
                                    int_norm[1] * t.first,
                                    int_norm[2] * t.first };

                  // Try again with this new edge
                  connect (nh, ico);
                  connect (ico.m.prev (ico.m.opposite (nh)), ico);
                  if (t.second != nullptr)
                    {
                      vcolor[ico.m.target (nh)].insert (t.second);
                      bfs_for_vert (ico, s, vcolor, evisited,
                                    voronoi_vertices);
                    }
                  // return;
                }
              else
                {
                  ico.m.point (ico.m.target (nh))
                      = K::Point_3{ int_norm[0] * t.first,
                                    int_norm[1] * t.first,
                                    int_norm[2] * t.first };
                  vcolor[ico.m.target (nh)].insert (curve_colors[0]);
                  vcolor[ico.m.target (nh)].insert (curve_colors[1]);
                  trace_out_edge (ico, curve_colors, nh, vcolor, evisited,
                                  voronoi_vertices);

                  trace_out_edge (ico, curve_colors, ico.m.opposite (h),
                                  vcolor, evisited, voronoi_vertices);
                  if (!CGAL::is_triangle_mesh (ico.m))
                    {
                      throw std::runtime_error (
                          "mesh is in invalid state after bfs");
                    }
                }
            }
        }
    }
}
