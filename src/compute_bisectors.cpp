#include "typedefs.h"
#include "icosphere.h"

inline Hyperbola const
compute_bisector_same (const Icosphere &ico, const Atom &atom, const Atom &it)
{

  // The coordinate system is scaled and translated such that the base atom is
  // at (0,0,0) and has radius 1. This cancels some of the terms in the full
  // equation for the quadric
  double x2 = (atom.x () - it.x ())
              / ico.radii_map (&atom);
  double y2 = (atom.y () - it.y ())
              / ico.radii_map (&atom);
  double z2 = (atom.z () - it.z ())
              / ico.radii_map (&atom);

  // Parameters for quadric. These are from
  //
  //  Zhongyin Hu, Xiang Li, Adarsh Krishnamurthy, Iddo Hanniel & Sara McMains
  //  (2017) Voronoi cells of non-general position spheres using the GPU,
  //  Computer-Aided Design and Applications, 14:5, 572-581,
  //  DOI: 10.1080/16864360.2016.1273576
  //
  // but I ran them through mathematica with r1=1, x1=0, y1=0, z1=0 to
  // to simplifiy them. There is also a typo in the paper; the definitions of E
  // and F were switched. Here they are corrected.
  double G = 2. * x2;
  double H = 2. * y2;
  double I = 2. * z2;
  double J = std::pow (x2, 2) + std::pow (y2, 2) + std::pow (z2, 2);

  return { 0., 0., 0., 0., 0., 0., G, H, I, J };
}

inline Hyperbola const
compute_bisector_diff (const Icosphere &ico, const Atom &atom, const Atom &it)
{

  // The coordinate system is scaled and translated such that the base atom is
  // at (0,0,0) and has radius 1. This cancels some of the terms in the full
  // equation for the quadric
  double x2 = (atom.x () - it.x ())
              / ico.radii_map (&atom);
  double y2 = (atom.y () - it.y ())
              / ico.radii_map (&atom);
  double z2 = (atom.z () - it.z ())
              / ico.radii_map (&atom);

  double r = (ico.radii_map (&atom)
              - ico.radii_map (&it))
             / ico.radii_map (&atom);

  // Helper variable for quadric
  double K = std::pow (x2, 2) + std::pow (y2, 2) + std::pow (z2, 2)
             - std::pow (r, 2);

  // Parameters for quadric. These are from
  //
  //  Zhongyin Hu, Xiang Li, Adarsh Krishnamurthy, Iddo Hanniel & Sara McMains
  //  (2017) Voronoi cells of non-general position spheres using the GPU,
  //  Computer-Aided Design and Applications, 14:5, 572-581,
  //  DOI: 10.1080/16864360.2016.1273576
  //
  // but I ran them through mathematica with r1=1, x1=0, y1=0, z1=0 to
  // to simplifiy them. There is also a typo in the paper; the definitions of E
  // and F were switched. Here they are corrected.
  double A = -4. * std::pow (x2, 2) + 4. * std::pow (r, 2);
  double B = -4. * std::pow (y2, 2) + 4. * std::pow (r, 2);
  double C = -4. * std::pow (z2, 2) + 4. * std::pow (r, 2);
  double D = -8. * x2 * y2;
  double E = -8. * x2 * z2;
  double F = -8. * y2 * z2;
  double G = -4. * x2 * K;
  double H = -4. * y2 * K;
  double I = -4. * z2 * K;
  double J = -std::pow (K, 2);

  //  printf("%s %i\n"
  //          "A=%f,B=%f,C=%f,D=%f,Ee=%f,F=%f,G=%f,H=%f,Ii=%f,J=%f\n",
  //          atom.element().c_str(), it.atom_serial_number(),
  //          A,B,C,D,E,F,G,H,I,J);

  return { A, B, C, D, E, F, G, H, I, J };
}

Hyperbola const
compute_bisector (const Icosphere &ico, const Atom &atom, const Atom &it)
{
  return ico.radii_map (&atom)
                 == ico.radii_map (&it)
             ? compute_bisector_same (ico, atom, it)
             : compute_bisector_diff (ico, atom, it);
}
