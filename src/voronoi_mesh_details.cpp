#include "awVd_ray.h"
#include "icosphere.h"
#include "ray_intersections.h"
#include "typedefs.h"
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Segment_3.h>
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <limits>
#include <utility>

#define ULP 100000

std::vector<Eigen::Matrix<NT, 3, 1> >
case_1 (const Eigen::Matrix<NT, 3, 5> &abcdf)
{
  NT R;

  // det(a,b,c)
  Eigen::Matrix<double, 3, 3> F_m;
  F_m << abcdf.col (0), abcdf.col (1), abcdf.col (2);
  const double F = F_m.determinant ();

  // det(f,b,c)
  Eigen::Matrix<double, 3, 3> F10_m;
  F10_m << abcdf.col (4), abcdf.col (1), abcdf.col (2);
  const double F10 = F10_m.determinant ();

  // det(-d,b,c)
  Eigen::Matrix<double, 3, 3> F11_m;
  F11_m << -1. * abcdf.col (3), abcdf.col (1), abcdf.col (2);
  const double F11 = F11_m.determinant ();

  // det(a,f,c)
  Eigen::Matrix<double, 3, 3> F20_m;
  F20_m << abcdf.col (0), abcdf.col (4), abcdf.col (2);
  const double F20 = F20_m.determinant ();

  // det(a,-d,c)
  Eigen::Matrix<double, 3, 3> F21_m;
  F21_m << abcdf.col (0), -1. * abcdf.col (3), abcdf.col (2);
  const double F21 = F21_m.determinant ();

  // det(a,b,f)
  Eigen::Matrix<double, 3, 3> F30_m;
  F30_m << abcdf.col (0), abcdf.col (1), abcdf.col (4);
  const double F30 = F30_m.determinant ();

  // det(a,b,-d)
  Eigen::Matrix<double, 3, 3> F31_m;
  F31_m << abcdf.col (0), abcdf.col (1), -1. * abcdf.col (3);
  const double F31 = F31_m.determinant ();

  // clang-format off

  const double alpha = std::pow(F11, 2) / std::pow(F, 2) +
                       std::pow(F21, 2) / std::pow(F, 2) +
                       std::pow(F31, 2) / std::pow(F, 2) - 1;

  const double beta = 2. * F10 * F11 / std::pow(F, 2) +
                      2. * F20 * F21 / std::pow(F, 2) +
                      2. * F30 * F31 / std::pow(F, 2) - 2.;

  const double gamma = std::pow(F10, 2) / std::pow(F, 2) +
                       std::pow(F20, 2) / std::pow(F, 2) +
                       std::pow(F30, 2) / std::pow(F, 2) - 1;
  // clang-format on

  double det = std::pow (beta, 2) - 4. * alpha * gamma;

  double x, y, z;

  // There is only one tangent sphere
  if (almost_equal (det, 0., 5))
    {
      R = -beta / (2. * alpha);
      x = F10 / F + R * F11 / F;
      y = F20 / F + R * F21 / F;
      z = F30 / F + R * F31 / F;
      return std::vector<Eigen::Matrix<NT, 3, 1> >{ { x, y, z } };
    }
  // There is no tangent sphere
  else if (det < 0)
    {
      return {};
    }
  // There are two tangent spheres
  else
    {
      double radius1 = (-beta - std::sqrt (det)) / (2. * alpha);
      double radius2 = (-beta + std::sqrt (det)) / (2. * alpha);
      // If there are 2 candidates, choose the one closest to the triangle
      double x1, y1, z1, x2, y2, z2;

      double Rmax = std::max (radius1, radius2);
      x1 = F10 / F + Rmax * F11 / F;
      y1 = F20 / F + Rmax * F21 / F;
      z1 = F30 / F + Rmax * F31 / F;

      double Rmin = std::min (radius1, radius2);
      x2 = F10 / F + Rmin * F11 / F;
      y2 = F20 / F + Rmin * F21 / F;
      z2 = F30 / F + Rmin * F31 / F;

      return std::vector<Eigen::Matrix<NT, 3, 1> >{ { x1, y1, z1 },
                                                    { x2, y2, z2 } };
    }
}

std::vector<Eigen::Matrix<NT, 3, 1> >
case_2 (Eigen::Matrix<NT, 3, 5> &abcdf)
{
  auto ColPivQR_abd
      = (abcdf.col (0), abcdf.col (1), abcdf.col (3)).colPivHouseholderQr
      ();
  int rank_abd = ColPivQR_abd.rank ();

  if (rank_abd == 3)
    {
      abcdf.col (2).swap (abcdf.col (3));
      return case_1 (abcdf);
    }

    auto ColPivQR_adc
        = (abcdf.col (0), abcdf.col (3), abcdf.col (2)).colPivHouseholderQr
        ();
    int rank_adc = ColPivQR_adc.rank ();

  if (rank_adc == 3)
    {
      abcdf.col (1).swap (abcdf.col (3));
      return case_1 (abcdf);
    }
  else
    {
      abcdf.col (0).swap (abcdf.col (3));
      return case_1 (abcdf);
    }
}

std::vector<Eigen::Matrix<NT, 3, 1> >
compute_voronoi_vertex_rigorous (const Icosphere &ico,
                                 const std::array<const Atom *, 3> atoms)
{

  // This solves a systems of equations to find a point equidistant from
  // the surface of 4 spheres. The coordinate sysem is scaled and
  // translated such that the base sphere is at (0,0,0) and has radius 1.
  NT x1 = ico.center ()[0];
  NT y1 = ico.center ()[1];
  NT z1 = ico.center ()[2];
  NT r1 = ico.radius ();

  // Scale and translate the coordinate system such that S1 has center (0,0,0)
  // and radius 1

  const Atom *atom2 = atoms[0];
  NT x2 = (atom2->x () - x1) / r1;
  NT y2 = (atom2->y () - y1) / r1;
  NT z2 = (atom2->z () - z1) / r1;

  NT r2 = G_atom_classifier.get_properties (*atom2).value () / r1;

  auto atom3 = atoms[1];
  NT x3 = (atom3->x () - x1) / r1;
  NT y3 = (atom3->y () - y1) / r1;
  NT z3 = (atom3->z () - z1) / r1;
  NT r3 = G_atom_classifier.get_properties (*atom3).value () / r1;

  auto atom4 = atoms[2];
  NT x4 = (atom4->x () - x1) / r1;
  NT y4 = (atom4->y () - y1) / r1;
  NT z4 = (atom4->z () - z1) / r1;
  NT r4 = G_atom_classifier.get_properties (*atom4).value () / r1;

  NT a1 = 2. * x2;
  NT b1 = 2. * y2;
  NT c1 = 2. * z2;
  NT d1 = 2. * r2 - 2.;
  NT f1 = 1. - std::pow (r2, 2) + std::pow (x2, 2) + std::pow (y2, 2)
          + std::pow (z2, 2);

  NT a2 = 2. * x3;
  NT b2 = 2. * y3;
  NT c2 = 2. * z3;
  NT d2 = 2. * r3 - 2.;
  NT f2 = 1. - std::pow (r3, 2) + std::pow (x3, 2) + std::pow (y3, 2)
          + std::pow (z3, 2);

  NT a3 = 2. * x4;
  NT b3 = 2. * y4;
  NT c3 = 2. * z4;
  NT d3 = 2. * r4 - 2.;
  NT f3 = 1. - std::pow (r4, 2) + std::pow (x4, 2) + std::pow (y4, 2)
          + std::pow (z4, 2);

  // clang-format off
  Eigen::Matrix<double, 3, 3> abc;
  abc << a1, b1, c1,
         a2, b2, c2,
         a3, b3, c3;

  Eigen::Matrix<double, 3, 4> abcd;
  abcd << a1, b1, c1, d1,
          a2, b2, c2, d2,
          a3, b3, c3, d3;

  Eigen::Matrix<double, 3, 5> abcdf;
  abcdf << a1, b1, c1, d1, f1,
          a2, b2, c2, d2, f2,
          a3, b3, c3, d3, f3;
  // clang-format on

  Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 3> > ColPivQR_abc
      = abc.colPivHouseholderQr ();
  int rank_abc = ColPivQR_abc.rank ();

  Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 4> > ColPivQR_abcd
      = abcd.colPivHouseholderQr ();
  int rank_abcd = ColPivQR_abcd.rank ();

  Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 5> > ColPivQR_abcdf
      = abcdf.colPivHouseholderQr ();
  int rank_abcdf = ColPivQR_abcdf.rank ();

  std::vector<Eigen::Matrix<NT, 3, 1> > ret;
  if (rank_abc == 3 && rank_abcd == 3 && rank_abcdf == 3)
    {
      //  printf ("case 1\n");

      auto v_verts = case_1 (abcdf);
      if (!v_verts.empty ())
        {
          double v1_to_s0 = std::sqrt (std::pow (v_verts[0][0], 2)
                                       + std::pow (v_verts[0][1], 2)
                                       + std::pow (v_verts[0][2], 2))
                            - 1.;
          double v1_to_s1 = std::sqrt (std::pow (v_verts[0][0] - x2, 2)
                                       + std::pow (v_verts[0][1] - y2, 2)
                                       + std::pow (v_verts[0][2] - z2, 2))
                            - r2;
          double v1_to_s2 = std::sqrt (std::pow (v_verts[0][0] - x3, 2)
                                       + std::pow (v_verts[0][1] - y3, 2)
                                       + std::pow (v_verts[0][2] - z3, 2))
                            - r3;
          double v1_to_s3 = std::sqrt (std::pow (v_verts[0][0] - x4, 2)
                                       + std::pow (v_verts[0][1] - y4, 2)
                                       + std::pow (v_verts[0][2] - z4, 2))
                            - r4;

          bool s0_s1 = almost_equal (v1_to_s0, v1_to_s1, ULP);
          bool s1_s2 = almost_equal (v1_to_s1, v1_to_s2, ULP);
          bool s2_s3 = almost_equal (v1_to_s2, v1_to_s3, ULP);
          bool s3_s0 = almost_equal (v1_to_s3, v1_to_s0, ULP);
          bool s2_s0 = almost_equal (v1_to_s2, v1_to_s0, ULP);
          bool s1_s3 = almost_equal (v1_to_s1, v1_to_s3, ULP);

          if (s0_s1 && s1_s2 && s2_s3 && s3_s0 && s2_s0 && s1_s3)
            {
              ret.push_back (v_verts[0]);
            }
          if (v_verts.size () > 1)
            {
              double v2_to_s0 = std::sqrt (std::pow (v_verts[1][0], 2)
                                           + std::pow (v_verts[1][1], 2)
                                           + std::pow (v_verts[1][2], 2))
                                - 1.;
              double v2_to_s1 = std::sqrt (std::pow (v_verts[1][0] - x2, 2)
                                           + std::pow (v_verts[1][1] - y2, 2)
                                           + std::pow (v_verts[1][2] - z2, 2))
                                - r2;
              double v2_to_s2 = std::sqrt (std::pow (v_verts[1][0] - x3, 2)
                                           + std::pow (v_verts[1][1] - y3, 2)
                                           + std::pow (v_verts[1][2] - z3, 2))
                                - r3;
              double v2_to_s3 = std::sqrt (std::pow (v_verts[1][0] - x4, 2)
                                           + std::pow (v_verts[1][1] - y4, 2)
                                           + std::pow (v_verts[1][2] - z4, 2))
                                - r4;

              s0_s1 = almost_equal (v2_to_s0, v2_to_s1, ULP);
              s1_s2 = almost_equal (v2_to_s1, v2_to_s2, ULP);
              s2_s3 = almost_equal (v2_to_s2, v2_to_s3, ULP);
              s3_s0 = almost_equal (v2_to_s3, v2_to_s0, ULP);
              s2_s0 = almost_equal (v2_to_s2, v2_to_s0, ULP);
              s1_s3 = almost_equal (v2_to_s1, v2_to_s3, ULP);


              bool v1_v2x = almost_equal(v_verts[0][0], v_verts[1][0], ULP);
              bool v1_v2y = almost_equal(v_verts[0][1], v_verts[1][1], ULP);
              bool v1_v2z = almost_equal(v_verts[0][2], v_verts[1][2], ULP);
              bool same_coords = v1_v2x && v1_v2y && v1_v2z;

              if (s0_s1 && s1_s2 && s2_s3 && s3_s0 && s2_s0 && s1_s3 && !same_coords)
                {
                  ret.push_back (v_verts[1]);
                }
            }
        }
    }
  else if (rank_abc == 2 && rank_abcd == 3 && rank_abcdf == 3)
    {
      printf ("case 2\n");
      return case_2 (abcdf);
    }
  else if (rank_abc == 2 && rank_abcd == 2 && rank_abcdf == 3)
    {
    }
  else if (rank_abc == 2 && rank_abcd == 2 && rank_abcdf == 2)
    {
    }
  else if (rank_abc == 1 && rank_abcd == 2 && rank_abcdf == 3)
    {
    }
  else if (rank_abc == 1 && rank_abcd == 2 && rank_abcdf == 2)
    {
    }
  else if (rank_abc == 1 && rank_abcd == 1 && rank_abcdf == 2)
    {
    }
  else if (rank_abc == 1 && rank_abcd == 1 && rank_abcdf == 1)
    {
    }
  return ret;
}

#define EPSILON 1e7

std::list<K::Point_3>
find_voronoi_curves_along_mesh_edge (
    const std::array<const Atom *, 2> &curve_colors, Icosphere &ico,
    const Mesh::Edge_index &edge)
{
  auto start = ico.m.vertex (edge, 0);
  auto _p = ico.m.point (start);
  Eigen::Matrix<NT, 3, 1> p;
  p << _p.x (), _p.y (), _p.z ();
  NT px = _p.x();
  NT py = _p.y();
  NT pz = _p.z();

  auto end = ico.m.vertex (edge, 1);
  auto _q = ico.m.point (end);
  Eigen::Matrix<NT, 3, 1> q;
  q << _q.x (), _q.y (), _q.z ();
  NT qx = _q.x();
  NT qy = _q.y();
  NT qz = _q.z();

  NT nx = py * qz - pz * qy;
  NT ny = pz * qx - px * qz;
  NT nz = px * qy - py * qx;

  NT x1 = ico.center ().x();
  NT y1 = ico.center ().y();
  NT z1 = ico.center ().z();
  NT r1 = ico.radius ();

  // Scale and translate the coordinate system such that S1 has center (0,0,0)
  // and radius 1

  std::list<K::Point_3> intersections;

  auto atom1 = curve_colors[0];
  auto atom2 = curve_colors[1];

  NT x2 = (atom1->x () - x1) / r1;
  NT y2 = (atom1->y () - y1) / r1;
  NT z2 = (atom1->z () - z1) / r1;
  NT r2 = G_atom_classifier.get_properties (*atom1).value () / r1;

  NT x3 = (atom2->x () - x1) / r1;
  NT y3 = (atom2->y () - y1) / r1;
  NT z3 = (atom2->z () - z1) / r1;
  NT r3 = G_atom_classifier.get_properties (*atom2).value () / r1;

  NT a1 = (NT)2. * x2;
  NT b1 = (NT)2. * y2;
  NT c1 = (NT)2. * z2;
  NT d1 = (NT)2. * r2 - (NT)2.;
  NT f1 = (NT)1. - (r2 * r2) + (x2 * x2) + (y2 * y2) + (z2 * z2);

  NT a2 = (NT)2. * x3;
  NT b2 = (NT)2. * y3;
  NT c2 = (NT)2. * z3;
  NT d2 = (NT)2. * r3 - (NT)2.;
  NT f2 = (NT)1. - (r3 * r3) + (x3 * x3) + (y3 * y3) + (z3 * z3);

  NT a3 = nx;
  NT b3 = ny;
  NT c3 = nz;
  NT d3 = 0.;
  NT f3 = 0.;

  NT R;

  // clang-format off
  NT F = a1 * b2 * c3 -
             a1 * b3 * c2 -
             a2 * b1 * c3 +
             a2 * b3 * c1 +
             a3 * b1 * c2 -
             a3 * b2 * c1;

  NT F10 = b1 * c2 * f3 -
               b1 * c3 * f2 -
               b2 * c1 * f3 +
               b2 * c3 * f1 +
               b3 * c1 * f2 -
               b3 * c2 * f1;

  NT F11 = -b1 * c2 * d3 +
                     b1 * c3 * d2 +
                     b2 * c1 * d3 -
                     b2 * c3 * d1 -
                     b3 * c1 * d2 +
                     b3 * c2 * d1;

  NT F20 =  -a1 * c2 * f3 +
                     a1 * c3 * f2 +
                     a2 * c1 * f3 -
                     a2 * c3 * f1 -
                     a3 * c1 * f2 +
                     a3 * c2 * f1;

  NT F21 = a1 * c2 * d3 -
               a1 * c3 * d2 -
               a2 * c1 * d3 +
               a2 * c3 * d1 +
               a3 * c1 * d2 -
               a3 * c2 * d1;

  NT F30 = a1 * b2 * f3 -
               a1 * b3 * f2 -
               a2 * b1 * f3 +
               a2 * b3 * f1 +
               a3 * b1 * f2 -
               a3 * b2 * f1;

  NT F31 =  -a1 * b2 * d3 +
                     a1 * b3 * d2 +
                     a2 * b1 * d3 -
                     a2 * b3 * d1 -
                     a3 * b1 * d2 +
                     a3 * b2 * d1;

  NT a = (F11 * F11) / (F * F) +
         (F21 * F21) / (F * F) +
         (F31 * F31) / (F * F) - (NT)1.;

  NT b = ((NT)2. * F10 * F11) / (F * F) +
         ((NT)2. * F20 * F21) / (F * F) +
         ((NT)2. * F30 * F31) / (F * F) - (NT)2.;

  NT c = (F10 * F10) / (F * F) +
         (F20 * F20) / (F * F) +
         (F30 * F30) / (F * F) - (NT)1.;
  // clang-format on

  NT det = (b * b) - (NT)4. * a * c;
  // There is only one tangent sphere
  if (almost_equal (det, 0., 500))
    {
      R = -b / (2. * a);
      NT x = F10 / F + R * F11 / F;
      NT y = F20 / F + R * F21 / F;
      NT z = F30 / F + R * F31 / F;
      //      const auto i = CGAL::intersection (
      //          K::Ray_3{ K::Point_3{ 0, 0, 0 }, K::Point_3{ x, y, z } },
      //          K::Segment_3{ ico.m.point (start), ico.m.point (end) });
      K::Plane_3 plane{ _q, _p, K::Point_3{ 0., 0., 0. } };

      K::Ray_2 ray{ plane.to_2d (K::Point_3{ 0., 0., 0. }),
                    plane.to_2d (K::Point_3{ x, y, z }) };
      K::Segment_2 seg{ plane.to_2d (_q), plane.to_2d (_p) };

      auto i = CGAL::intersection (ray, seg);
      NT v_to_s0
          = std::sqrt (std::pow (x, 2) + std::pow (z, 2) + std::pow (y, 2))
            - 1.;
      NT v_to_s1 = std::sqrt (std::pow (x - x2, 2) + std::pow (z - z2, 2)
                              + std::pow (y - y2, 2))
                   - r2;
      NT v_to_s2 = std::sqrt (std::pow (x - x3, 2) + std::pow (z - z3, 2)
                              + std::pow (y - y3, 2))
                   - r3;

      bool s0_s1 = almost_equal (v_to_s0, v_to_s1, ULP);
      bool s1_s2 = almost_equal (v_to_s1, v_to_s2, ULP);
      bool s2_s0 = almost_equal (v_to_s2, v_to_s0, ULP);
      if (i && s0_s1 && s1_s2 && s2_s0)
        {
          intersections.push_back (K::Point_3{ x, y, z });
        }
    }
  // There are two tangent spheres
  else
    {
      NT radius1 = (-b - std::sqrt (det)) / (2. * a);
      NT radius2 = (-b + std::sqrt (det)) / (2. * a);

      NT x_1 = F10 / F + radius1 * F11 / F;
      NT y_1 = F20 / F + radius1 * F21 / F;
      NT z_1 = F30 / F + radius1 * F31 / F;

      NT x_2 = F10 / F + radius2 * F11 / F;
      NT y_2 = F20 / F + radius2 * F21 / F;
      NT z_2 = F30 / F + radius2 * F31 / F;

      K::Plane_3 plane{ _q, _p, K::Point_3{ 0., 0., 0. } };

      K::Ray_2 ray1{ plane.to_2d (K::Point_3{ 0., 0., 0. }),
                     plane.to_2d (K::Point_3{ x_1, y_1, z_1 }) };
      K::Ray_2 ray2{ plane.to_2d (K::Point_3{ 0., 0., 0. }),
                     plane.to_2d (K::Point_3{ x_2, y_2, z_2 }) };
      K::Segment_2 seg{ plane.to_2d (_q), plane.to_2d (_p) };

      auto i1 = CGAL::intersection (ray1, seg);
      NT v1_to_s0 = std::sqrt (std::pow (x_1, 2) + std::pow (z_1, 2)
                               + std::pow (y_1, 2))
                    - 1.;
      NT v1_to_s1 = std::sqrt (std::pow (x_1 - x2, 2) + std::pow (z_1 - z2, 2)
                               + std::pow (y_1 - y2, 2))
                    - r2;
      NT v1_to_s2 = std::sqrt (std::pow (x_1 - x3, 2) + std::pow (z_1 - z3, 2)
                               + std::pow (y_1 - y3, 2))
                    - r3;

      bool s0_s1 = almost_equal (v1_to_s0, v1_to_s1, ULP);
      bool s1_s2 = almost_equal (v1_to_s1, v1_to_s2, ULP);
      bool s2_s0 = almost_equal (v1_to_s2, v1_to_s0, ULP);

      if (i1 && s0_s1 && s1_s2 && s2_s0)
        {
          intersections.push_back (K::Point_3{ x_1, y_1, z_1 });
        }
      auto i2 = CGAL::intersection (ray2, seg);
      NT v2_to_s0 = std::sqrt (std::pow (x_2, 2) + std::pow (z_2, 2)
                               + std::pow (y_2, 2))
                    - 1.;
      NT v2_to_s1 = std::sqrt (std::pow (x_2 - x2, 2) + std::pow (z_2 - z2, 2)
                               + std::pow (y_2 - y2, 2))
                    - r2;
      NT v2_to_s2 = std::sqrt (std::pow (x_2 - x3, 2) + std::pow (z_2 - z3, 2)
                               + std::pow (y_2 - y3, 2))
                    - r3;

      s0_s1 = almost_equal (v2_to_s0, v2_to_s1, ULP);
      s1_s2 = almost_equal (v2_to_s1, v2_to_s2, ULP);
      s2_s0 = almost_equal (v2_to_s2, v2_to_s0, ULP);

      bool x1_x2 = almost_equal (x_1, x_2, ULP);
      bool y1_y2 = almost_equal (y_1, y_2, ULP);
      bool z1_z2 = almost_equal (z_1, z_2, ULP);
      bool same_coords = x1_x2 && y1_y2 && z1_z2;

      if (i2 && s0_s1 && s1_s2 && s2_s0 && !same_coords)
        {
          intersections.emplace_front (K::Point_3{ x_2, y_2, z_2 });
        }
    }

  // assert (intersections.size () < 2);
  return intersections;
}
