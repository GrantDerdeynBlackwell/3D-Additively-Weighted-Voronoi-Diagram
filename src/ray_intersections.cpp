#include "icosphere.h"
#include "typedefs.h"
#include <CGAL/number_utils.h>
#include <limits>

bool
almost_equal (double x, double y, std::size_t ulp)
{
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::fabs (x - y) <= std::numeric_limits<double>::epsilon ()
                                  * std::fabs (x + y) * ulp
         // unless the result is subnormal
         || std::fabs (x - y) < std::numeric_limits<double>::min ();
}

double
compute_scalar_dist_different_sizes (
    const Ray &ray, const std::pair<const Atom *const, Hyperbola> &hyperbola,
    const bool bigger)
{
  NT t = G_double_to_NT (std::numeric_limits<double>::max ());

  // The icosphere has radius 1 already, so these don't need to be normalized
  const NT nx = ray[0];
  const NT ny = ray[1];
  const NT nz = ray[2];

  // Compute the ray distance from the origin. This makes the calculations with
  // overlapping spheres easier to handle.
  const NT Ox = 0;
  const NT Oy = 0;
  const NT Oz = 0;

  const NT A = hyperbola.second[0];
  const NT B = hyperbola.second[1];
  const NT C = hyperbola.second[2];
  const NT D = hyperbola.second[3];
  const NT E = hyperbola.second[4];
  const NT F = hyperbola.second[5];
  const NT G = hyperbola.second[6];
  const NT H = hyperbola.second[7];
  const NT I = hyperbola.second[8];
  const NT J = hyperbola.second[9];

  // I am confident in these. I calculated them myself and confirmed them
  // with Practical Ray Tracing in C by Lindley, Craig pp. 86-89.

  // clang-format off
  const NT a = A * CGAL::square (nx) + D * nx * ny +
                   B * CGAL::square(ny) + E * nx * nz +
                   C * CGAL::square (nz) + F * ny * nz;

  const NT b = 2 * A * nx * Ox + D * ny * Ox + D * nx * Oy + G * nx +
                   2 * B * ny * Oy + E * nz * Ox + E * nx * Oz + H * ny +
                   2 * C * nz * Oz + F * nz * Oy + F * ny * Oz + I * nz;

  const NT c = A * CGAL::square (Ox) + D * Ox * Oy + G * Ox +
                   B * CGAL::square (Oy) + E * Ox * Oz + H * Oy +
                   C * CGAL::square (Oz) + F * Oy * Oz + I * Oz + J;
  // clang-format on

  const NT det = CGAL::square (b) - 4. * a * c;

  // Only change things if there is an intersection
  if (almost_equal (det, 0., 100))
    {
      NT t1 = -1. * b / (2. * a);
      return t1 > 0 ? t1 : t;
      // only keep this sheet if it is the closest encountered, but not
      // negative

      // only update the color if the distance was changed
    }
  else if (det >= 0.)
    {
      const NT t1 = (-1. * b - CGAL::sqrt (det)) / (2. * a);
      const NT t2 = (-1. * b + CGAL::sqrt (det)) / (2. * a);

      // If the other atom is bigger than the reference atom, choose the
      // more negative bisector. Discard this bisector if it has a negative
      // scalar distance.
      if (bigger)
        {
          if (t1 > 0. && t2 > 0.)
            {
              t = CGAL::min (t1, t2);
            }
          else if (t2 > 0.)
            {
              t = t2;
            }
        }
      // If the reference atom is bigger than the other atom, choose the
      // more postive bisector. Discard this bisector if it has a negative
      // scalar distance.
      else
        {
          if (t1 > 0. && t2 > 0.)
            {
              t = CGAL::max (t1, t2);
            }
          else if (t1 > 0.)
            {
              t = t1;
            }
        }
    }
  return t;
}

NT
compute_scalar_dist_same_sizes (
    const Ray &ray, const std::pair<const Atom *const, Hyperbola> &hyperbola)
{

  // The icosphere has radius 1 already, so these don't need to be normalized
  const NT nx = ray[0];
  const NT ny = ray[1];
  const NT nz = ray[2];

  // Compute the ray distance from the origin. This makes the calculations with
  // overlapping spheres easier to handle.

  const NT G = hyperbola.second[6];
  const NT H = hyperbola.second[7];
  const NT I = hyperbola.second[8];
  const NT J = hyperbola.second[9];

  NT t = -1. * J / (nx * G + ny * H + nz * I);

  return t > 0. ? t : G_double_to_NT (std::numeric_limits<double>::max ());
}

std::pair<NT, const Atom *>
compute_scalar_dist (const Ray &ray, const Icosphere &ico)
{
  const auto &atom = ico.atom ();
  const auto &h = ico.h ();
  NT t = G_double_to_NT (std::numeric_limits<double>::max ());
  NT d_curr = ico.cutoff () + ico.radii_map (&atom);
  const Atom *h_curr = nullptr;

  for (const auto &hyperbola : h)
    {
      if (ico.radii_map (&atom) == ico.radii_map (hyperbola.first))
        {
          t = compute_scalar_dist_same_sizes (ray, hyperbola);
        }
      else
        {
          bool bigger
              = ico.radii_map (&atom) < ico.radii_map (hyperbola.first);
          t = compute_scalar_dist_different_sizes (ray, hyperbola, bigger);
        }
      // only keep this sheet if it is the closest yet encountered.
      d_curr = std::min (d_curr, t);

      // only update the color if the distance was changed.
      h_curr = d_curr == t ? hyperbola.first : h_curr;
    }
  return std::make_pair (d_curr, h_curr);
}
