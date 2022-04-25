#ifndef G_ICOSPHERE
#define G_ICOSPHERE
#include "typedefs.h"
#include "parser.h"
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <cmath>
#include <cstddef>

class Icosphere
{
  const K::Point_3 _center;
  const NT _radius;
  const std::size_t _subdivisions;
  const Hyperbola_map &_h;
  const Atom &_atom;

public:
  Mesh m;

  Icosphere (const Atom &atom, const std::size_t subdivisions,
             const Hyperbola_map &h)
      : _center{ atom.x (), atom.y (), atom.z () },
        _radius (G_atom_classifier.get_properties (atom).value ()),
        _subdivisions (subdivisions), _h (h), _atom (atom)
  {
    m.reserve (20 * std::pow (4, subdivisions),
               40 * std::pow (4, subdivisions),
               30 * std::pow (4, subdivisions));

    // create a property map to store the 'color' of each vertex of the mesh
    Mesh::Property_map<Mesh::Edge_index, bool> evisited;
    bool ecreated;
    boost::tie (evisited, ecreated)
        = m.add_property_map<Mesh::Edge_index, bool> ("e:visited", false);
    assert (ecreated);

    Mesh::Property_map<vertex_descriptor, std::set<const Atom *> > vcolor;
    bool vcreated;
    boost::tie (vcolor, vcreated)
        = m.add_property_map<vertex_descriptor, std::set<const Atom *> > (
            "v:atom", std::set<const Atom *> ());
    assert (vcreated);

    CGAL::make_icosahedron (m, K::Point_3 (0., 0., 0.));


    // clang-format off
    CGAL::Aff_transformation_3<K> rotx
    (
      1.,             0.,              0.,
      0., std::cos (0.1), -std::sin (0.1),
      0., std::sin (0.1),  std::cos (0.1),
      1.
    );

    CGAL::Aff_transformation_3<K> roty
    (
      std::cos(0.1),  0.,   std::sin(0.1),
      0.,             1.,              0.,
      -std::sin(0.1), 0.,  std::cos (0.1),
      1.
    );

    CGAL::Aff_transformation_3<K> rotz
    (
      std::cos (0.1), -std::sin (0.1), 0.,
      std::sin (0.1),  std::cos (0.1), 0.,
      0.,                          0., 1.,
      1.
    );
    // clang-format on

    // The mesh starts align to each axis. We give it a little nudge to prevent
    // annoying numerical errors.
    CGAL::Polygon_mesh_processing::transform(rotx,m);
    CGAL::Polygon_mesh_processing::transform(roty,m);
    CGAL::Polygon_mesh_processing::transform(rotz,m);

    // Add the icosahedron to the mesh, coloring each vertex null for now
    // subdivide the mesh

    if (subdivisions > 0)
      {
        CGAL::Subdivision_method_3::Loop_subdivision (
            m, CGAL::parameters::number_of_iterations (subdivisions));
      }
  }

  const K::Point_3
  center () const
  {
    return this->_center;
  }

  double
  radius () const
  {
    return this->_radius;
  }

  const Hyperbola_map &
  h () const
  {
    return this->_h;
  }

  const Atom &
  atom () const
  {
    return this->_atom;
  }
  std::size_t
  subdivisions () const
  {
    return this->_subdivisions;
  }
};

#endif
