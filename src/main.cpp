#include "awVd_ray.h"
#include "power.h"
#include "typedefs.h"
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <ESBTL/PDB.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/default.h>
#include <ESBTL/occupancy_handlers.h>
#include <Eigen/Eigen>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <string>
namespace ESBTL
{
std::ostream &
operator<< (std::ostream &os, const Atom &atm)
{
  os << ESBTL::PDB::get_atom_pdb_format (atm);
  return os;
}
}

namespace po = boost::program_options;

int
main (int argc, const char **argv)
{
  po::options_description desc ("Allowed options");
  desc.add_options () ("help", "produce help message") (
      "input", po::value<std::string> ()->required (),
      "pdb file to analyze") ("subdivisions", po::value<int> ()->required (),
                              "number of subdivisions") (
      "vp", po::value<std::string> ()->implicit_value (""),
      "write power volumes to pdb") (
      "vv", po::value<std::string> ()->implicit_value (""),
      "write awVd volumes to pdb") ("vd", po::value<std::string> (),
                                    "output % difference in volumes to pdb") (
      "op", po::value<std::string> ()->implicit_value (""),
      "write the power cell-sphere overlaps to pdb") (
      "ov", po::value<std::string> ()->implicit_value (""),
      "write the awVd cell-sphere overlaps to pdb") (
      "od", po::value<std::string> ()->implicit_value (""),
      "write the % difference in cell-sphere overlaps to pdb") (
      "ap", po::value<std::string> ()->implicit_value (""),
      "output power diagram cell surface areas to pdb") (
      "av", po::value<std::string> ()->implicit_value (""),
      "output awVd cell surface areas to pdb") (
      "ad", po::value<std::string> ()->implicit_value (""),
      "output % difference in cell surface areas to pdb") (
      "ip", po::value<std::string> ()->implicit_value (""),
      "output the interfacial power cell surface areas to pdb") (
      "iv", po::value<std::string> ()->implicit_value (""),
      "output the interfacial awVd cell surface areas to pdb") (
      "id", po::value<std::string> ()->implicit_value (""),
      "output % difference in interfacial cell surface areas to pdb") (
      "k", po::value<std::string> ()->implicit_value (""),
      "output maximum Gaussian curvature of each cell to pdb") (
      "c", po::value<std::string> (),
      "name of csv") ("mp", po::value<bool> (),
                      "write a COFF mesh for each power diagram cell") (
      "mv", po::value<bool> (), "write a COFF mesh for each awVd cell") (
      "ra", po::value<std::vector<std::string> > ()->multitoken (),
      "include residues") (
      "rs", po::value<std::vector<std::string> > ()->multitoken (),
      "exclude residues") ("verts,v", po::bool_switch (),
                           "output the voronoi vertices to stdout");

  po::positional_options_description p;
  p.add ("input", 1);
  p.add ("subdivisions", 1);

  po::variables_map vm;
  po::store (po::command_line_parser (argc, argv)
                 .options (desc)
                 .positional (p)
                 .run (),
             vm);
  po::notify (vm);

  ESBTL::PDB_line_selector sel;

  std::vector<System> systems;
  ESBTL::All_atom_system_builder<System> builder (systems,
                                                  sel.max_nb_systems ());

  std::string filename = vm["input"].as<std::string> ();
  const int density = vm["subdivisions"].as<int> ();

  std::vector<std::string> default_residues{
    "DA",  "DG",  "DT",  "DC",  "ALA", "ARG", "ASN", "ASP", "ASX",
    "CYS", "CYX", "GLU", "GLN", "GLX", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
  };

  std::vector<std::string> include;
  std::vector<std::string> exclude;

  if (vm.count ("ra"))
    {
      include = vm["ra"].as<std::vector<std::string> > ();
    }
  if (vm.count ("rs"))
    {
      exclude = vm["rs"].as<std::vector<std::string> > ();
    }

  std::sort (include.begin (), include.end ());
  std::sort (exclude.begin (), exclude.end ());
  std::sort (default_residues.begin (), default_residues.end ());

  std::set<std::string> residues;

  std::set_difference (default_residues.begin (), default_residues.end (),
                       exclude.begin (), exclude.end (),
                       std::inserter (residues, residues.end ()));

  residues.insert (include.begin (), include.end ());

  if (ESBTL::read_a_pdb_file (filename, sel, builder,
                              Accept_none_occupancy_policy ()))
    {
      if (systems.empty () || systems[0].has_no_model ())
        {
          std::cerr << "No atoms found" << std::endl;
          return EXIT_FAILURE;
        }

      Model &model = *systems[0].models_begin ();

      std::map<std::string, std::ofstream> outs;

      std::ofstream csv;
      if (vm.count ("c"))
        {
          csv = std::ofstream (vm["c"].as<std::string> ());
        }
      else
        {
          csv = std::ofstream ("volumes.csv");
        }

      if (vm.count ("vp"))
        {
          outs["vp"] = std::ofstream (vm["vp"].as<std::string> () == ""
                                          ? "power_volumes.pdb"
                                          : vm["vp"].as<std::string> ());
        }
      if (vm.count ("vv"))
        {
          outs["vv"] = std::ofstream (vm["vv"].as<std::string> () == ""
                                          ? "awVd_volumes.pdb"
                                          : vm["vv"].as<std::string> ());
        }
      if (vm.count ("vd"))
        {
          outs["vd"] = std::ofstream (vm["vd"].as<std::string> () == ""
                                          ? "diff_volumes.pdb"
                                          : vm["vd"].as<std::string> ());
        }
      if (vm.count ("op"))
        {
          outs["op"] = std::ofstream (vm["op"].as<std::string> () == ""
                                          ? "power_overlap_volumes.pdb"
                                          : vm["op"].as<std::string> ());
        }
      if (vm.count ("ov"))
        {
          outs["ov"] = std::ofstream (vm["ov"].as<std::string> () == ""
                                          ? "awVd_overlap_volumes.pdb"
                                          : vm["ov"].as<std::string> ());
        }
      if (vm.count ("od"))
        {
          outs["od"] = std::ofstream (vm["od"].as<std::string> () == ""
                                          ? "diff_overlap_volumes.pdb"
                                          : vm["od"].as<std::string> ());
        }
      if (vm.count ("ap"))
        {
          outs["ap"] = std::ofstream (vm["ap"].as<std::string> () == ""
                                          ? "power_surface_areas.pdb"
                                          : vm["ap"].as<std::string> ());
        }
      if (vm.count ("av"))
        {
          outs["av"] = std::ofstream (vm["av"].as<std::string> () == ""
                                          ? "awVd_surface_areas.pdb"
                                          : vm["av"].as<std::string> ());
        }
      if (vm.count ("ad"))
        {
          outs["ad"] = std::ofstream (vm["ad"].as<std::string> () == ""
                                          ? "diff_surface_areas.pdb"
                                          : vm["ad"].as<std::string> ());
        }
      if (vm.count ("ip"))
        {
          outs["ip"] = std::ofstream (vm["ip"].as<std::string> () == ""
                                          ? "power_interface_areas.pdb"
                                          : vm["ip"].as<std::string> ());
        }
      if (vm.count ("iv"))
        {
          outs["iv"] = std::ofstream (vm["iv"].as<std::string> () == ""
                                          ? "awVd_interface_areas.pdb"
                                          : vm["iv"].as<std::string> ());
        }
      if (vm.count ("id"))
        {
          outs["id"] = std::ofstream (vm["id"].as<std::string> () == ""
                                          ? "diff_interface_areas.pdb"
                                          : vm["id"].as<std::string> ());
        }
      if (vm.count ("k"))
        {
          outs["k"] = std::ofstream (vm["k"].as<std::string> () == ""
                                         ? "max_curvature.pdb"
                                         : vm["k"].as<std::string> ());
        }
      std::map<const Atom *, std::map<std::string, double> > prop_map;

      // A map associating a quadruple of atoms with an empty sphere tangent to
      // all 4 atoms. Key: {cx, cy, cz, r}
      Voronoi_map voronoi_vertices;

      for (Atom_Iter it = model.atoms_begin (); it != model.atoms_end (); ++it)
        {
          std::array<double, 5> awVd_vol{ 0., 0., 0., 0., 0. };
          if (residues.find (it->residue_name ()) != residues.end ())
            {
              awVd_vol = find_neighbors (model, *(it), density, vm, residues,
                                         voronoi_vertices);
            }
          prop_map[&*it]["vv"] = awVd_vol[0];
          prop_map[&*it]["ov"] = awVd_vol[1];
          prop_map[&*it]["av"] = awVd_vol[2];
          prop_map[&*it]["iv"] = awVd_vol[3];
          prop_map[&*it]["k"] = awVd_vol[4];
        }

      Rt T;

      power (model, T);

      for (Rt::Vertex_handle vh : T.finite_vertex_handles ())
        {
          std::array<double, 4> power_vol{ 0., 0., 0., 0. };
          if (residues.find (vh->info ()->residue_name ()) != residues.end ())
            {
              power_vol = subdivide (vh, T, vm, residues);
            }
          prop_map[vh->info ()]["vp"] = power_vol[0];
          prop_map[vh->info ()]["op"] = power_vol[1];
          prop_map[vh->info ()]["ap"] = power_vol[2];
          prop_map[vh->info ()]["ip"] = power_vol[3];
        }

      outs["c"]
          << "atom,awVd_volume,awVd_overlap_volume,awVd_surface_area,awVd_"
             "interfacial_surface_area,maximum_gaussian_curvature,power_"
             "volume,power_overlap_volume,power_surface_area,power_"
             "interfacial_surface_area,%diff_volume,%diff_overlap_"
             "volume,%diff_surface_area,%diff_interfacial_surface_area"
          << std::endl;
      for (auto &atom_prop : prop_map)
        {
          auto &prop = atom_prop.second;
          auto atom = *atom_prop.first;
          prop["vd"] = prop["vp"] != 0.
                           ? ((prop["vv"] - prop["vp"]) / prop["vd"]) * 100.
                           : 0.;
          prop["od"] = prop["op"] != 0.
                           ? ((prop["ov"] - prop["op"]) / prop["op"]) * 100.
                           : 0.;
          prop["ad"] = prop["ap"] != 0.
                           ? ((prop["av"] - prop["ap"]) / prop["ap"]) * 100.
                           : 0.;
          prop["id"] = prop["ip"] != 0.
                           ? ((prop["iv"] - prop["ip"]) / prop["ip"]) * 100.
                           : 0.;
          csv << atom.atom_name () << "_" << atom.residue_name () << "_"
              << atom.atom_serial_number ();
          for (auto &elem : outs)
            {
              const auto &key = elem.first;
              atom.temperature_factor () = prop[key];
              if (outs[key].is_open ())
                {
                  outs[key] << atom << std::endl;
                }
              csv << "," << prop[key];
            }
        }
      if (vm["verts"].as<bool> ())
        {
          for (const auto &elem : voronoi_vertices)
            {
              const auto &atoms = elem.first;
              const auto &point = elem.second;
              for (const auto &atom : atoms)
                {
                  std::cout << atom->atom_serial_number () << ",";
                }
              std::cout << "\t" << point[0] << "," << point[1] << ","
                        << point[2] << std::endl;
            }
        }
    }
}
