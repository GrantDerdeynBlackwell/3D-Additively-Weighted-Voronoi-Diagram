#include "awVd_ray.h"
#include "parser.h"
#include "power.h"
#include "typedefs.h"
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <ESBTL/PDB.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/default.h>
#include <ESBTL/occupancy_handlers.h>
#include <Eigen/Eigen>
#include <cstdlib>
#include <filesystem>
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

void
parser (int argc, const char **argv, Cmd_line_options &options)
{
  for (int i = 3; i < argc; i++)
    {
      if (argv[i][0] == '-')
        {
          switch (argv[i][1])
            {
            case 'v':
              // option -v* is defined
              switch (argv[i][2])
                {
                case 'p':
                  // option -vp is defined. Output power diagram volumes to pdb
                  options.output_power = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-')
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto fname = arg.filename ();
                          auto parent = arg.parent_path ();
                          options.power_pdb_fname
                              = arg.replace_extension (".pdb");
                          options.power_pdb_overlap_fname
                              = arg.replace_filename (
                                  arg.stem ()
                                      .append ("_overlap")
                                      .replace_extension (".pdb"));
                        }
                    }
                  break;
                case 'v':
                  // option -vv is defined. Output awVd volumes to pdb
                  options.output_awVd = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-' && i < argc)
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto fname = arg.filename ();
                          auto parent = arg.parent_path ();
                          options.awVd_pdb_fname
                              = arg.replace_extension (".pdb");
                          options.awVd_pdb_overlap_fname
                              = arg.replace_filename (
                                  arg.stem ()
                                      .append ("_overlap")
                                      .replace_extension (".pdb"));
                        }
                    }
                  break;
                case 'd':
                  // option -vd is defined. Output the difference in the awVd
                  // and power volumes to pdb
                  options.output_diff = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-' && i < argc)
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto fname = arg.filename ();
                          auto parent = arg.parent_path ();
                          options.diff_pdb_fname
                              = arg.replace_extension (".pdb");
                          options.diff_pdb_overlap_fname
                              = arg.replace_filename (
                                  arg.stem ()
                                      .append ("_overlap")
                                      .replace_extension (".pdb"));
                        }
                    }
                  break;
                }
              break;
            case 'a':
              // option -a* is defined
              switch (argv[i][2])
                {
                case 'p':
                  // option -ap is defined. Output power diagram surface areas
                  // to pdb
                  options.output_power_areas = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-')
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto fname = arg.filename ();
                          auto parent = arg.parent_path ();
                          options.pow_area_pdb_fname
                              = arg.replace_extension (".pdb");
                          options.pow_interface_area_fname
                              = arg.replace_filename (
                                  arg.stem ()
                                      .append ("_interface")
                                      .replace_extension (".pdb"));
                        }
                    }
                  break;
                case 'v':
                  // option -av is defined. Output awVd surface areas to pdb
                  options.output_awVd_areas = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-' && i < argc)
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto fname = arg.filename ();
                          auto parent = arg.parent_path ();
                          options.awVd_area_pdb_fname
                              = arg.replace_extension (".pdb");
                          options.awVd_interface_area_fname
                              = arg.replace_filename (
                                  arg.stem ()
                                      .append ("_interface")
                                      .replace_extension (".pdb"));
                        }
                    }
                  break;
                case 'd':
                  // option -ad is defined. Output the difference in the awVd
                  // and power volumes to pdb
                  options.output_areas_diff = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-' && i < argc)
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto fname = arg.filename ();
                          auto parent = arg.parent_path ();
                          options.diff_area_pdb_fname
                              = arg.replace_extension (".pdb");
                          options.diff_area_pdb_interface_fname
                              = arg.replace_filename (
                                  arg.stem ()
                                      .append ("_interace")
                                      .replace_extension (".pdb"));
                        }
                    }
                  break;
                }
              break;
            case 'm':
              // option -m* is defined
              switch (argv[i][2])
                {
                case 'p':
                  // option -mp is defined. Output power diagram meshes
                  options.mesh_power = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-' && i < argc)
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto parent = arg.parent_path ();
                          options.power_mesh_dir = parent;
                        }
                    }
                  break;
                case 'v':
                  // option -mv is defined. Output awVd meshes
                  options.mesh_awVd = true;
                  if (i < argc - 1)
                    {
                      if (argv[i + 1][0] != '-' && i < argc)
                        {
                          // a filename has been given. Override the default
                          // value
                          std::filesystem::path arg
                              = (std::filesystem::path)argv[i + 1];
                          auto parent = arg.parent_path ();
                          options.awVd_mesh_dir = parent;
                        }
                    }
                  break;
                }
              break;
            case 'k':
              // option -k is degined. Output maximum curvature to pdb
              options.output_curvature = true;
              if (i < argc - 1)
                {
                  if (argv[i + 1][0] != '-' && i < argc)
                    {
                      // a filename has been given. Override the default value
                      std::filesystem::path arg
                          = (std::filesystem::path)argv[i + 1];
                      auto parent = arg.parent_path ();
                      options.curvature_pdb_fname = arg;
                    }
                }
              break;
            case 'c':
              // Option -c is defined. Output all the parameters to a csv
              options.output_csv = true;
              if (i < argc - 1)
                {
                  if (argv[i + 1][0] != '-' && i < argc)
                    {
                      // a filename has been given. Override the default value
                      std::filesystem::path arg
                          = (std::filesystem::path)argv[i + 1];
                      auto parent = arg.parent_path ();
                      options.csv = arg;
                    }
                }
              break;
            }
        }
    }
}

int
main (int argc, const char **argv)
{
  ESBTL::PDB_line_selector sel;

  std::vector<System> systems;
  ESBTL::All_atom_system_builder<System> builder (systems,
                                                  sel.max_nb_systems ());

  std::string mesh_filename;
  Cmd_line_options options;
  if (argc < 2)
    {
      std::cerr << "Please provide a filename" << std::endl;
      return EXIT_FAILURE;
    }
  if (argc < 3)
    {
      if ((std::string)argv[1] == "-h")
        {
          // Option -h is defined. Print help text
          printf (
              "Usage: vdos FILENAME SUBDIVISIONS [OPTIONS]\n"
              "Compute the Voronoi diagram of spheres for a pdb\n"
              "Example: vdos example.pdb 2\n"
              "Output control:\n"
              "-vp   output power diagram volumes to pdb\n"
              "-vv   output awVd volumes to pdb\n"
              "-vd   output %% difference in volumes to pdb\n"
              "-ap   output power diagram surface areas to pdb\n"
              "-av   output awVd surface areas to pdb\n"
              "-ad   output %% difference in surface areas to pdb\n"
              "-k    output maximum Gaussian curvature of each cell to pdb\n"
              "-c    output all of the above to a .csv\n"
              "-mp   write a COFF mesh for each power diagram cell\n"
              "-mv   write a COFF mesh for each awVd cell\n"
              "-h    print this help text\n");
          return EXIT_FAILURE;
        }
      std::cerr << "Please give the number of subdivisions" << std::endl;
      return EXIT_FAILURE;
    }

  std::string filename = argv[1];

  const int density = atoi (argv[2]);
  parser (argc, argv, options);

  if (ESBTL::read_a_pdb_file (filename, sel, builder,
                              Accept_none_occupancy_policy ()))
    {
      if (systems.empty () || systems[0].has_no_model ())
        {
          std::cerr << "No atoms found" << std::endl;
          return EXIT_FAILURE;
        }

      Model &model = *systems[0].models_begin ();

      std::array<std::ofstream, 13> outs;

      auto csv = std::ofstream (options.csv.string ());
      if (options.output_awVd)
        {
          printf ("writing awVd volumes to pdb\n");
          outs[0] = std::ofstream (options.awVd_pdb_fname.string ());
          outs[1] = std::ofstream (options.awVd_pdb_overlap_fname.string ());
        }
      if (options.output_awVd_areas)
        {
          printf ("writing awVd surface areas to pdb\n");
          outs[2] = std::ofstream (options.awVd_area_pdb_fname.string ());
          outs[3]
              = std::ofstream (options.awVd_interface_area_fname.string ());
        }
      if (options.output_curvature)
        {
          printf ("writing maximum_curvature to pdb\n");
          outs[4] = std::ofstream (options.curvature_pdb_fname.string ());
        }
      if (options.output_power)
        {
          printf ("writing power volumes to pdb\n");
          outs[5] = std::ofstream (options.power_pdb_fname.string ());
          outs[6] = std::ofstream (options.power_pdb_overlap_fname.string ());
        }
      if (options.output_power_areas)
        {
          printf ("writing power surface areas to pdb\n");
          outs[7] = std::ofstream (options.pow_area_pdb_fname.string ());
          outs[8] = std::ofstream (options.pow_interface_area_fname.string ());
        }
      if (options.output_diff)
        {
          printf ("writing difference in volumes to pdb\n");
          outs[9] = std::ofstream (options.diff_pdb_fname.string ());
          outs[10] = std::ofstream (options.diff_pdb_overlap_fname.string ());
        }
      if (options.output_areas_diff)
        {
          printf ("writing difference in volumes to pdb\n");
          outs[11] = std::ofstream (options.diff_area_pdb_fname.string ());
          outs[12] = std::ofstream (
              options.diff_area_pdb_interface_fname.string ());
        }

      // Keep track of the properties of each atom.
      //
      // awVd
      // ---------------------------------------------
      // prop_map[atom][0] = volume
      // prop_map[atom][1] = overlap_volume
      // prop_map[atom][2] = surface_area
      // prop_map[atom][3] = interfacial_surface_area
      // prop_map[atom][4] = maximum_curvature
      //
      // power
      // ---------------------------------------------
      // prop_map[atom][5] = volume
      // prop_map[atom][6] = overlap_volume
      // prop_map[atom][7] = surface_area
      // prop_map[atom][8] = interfacial_surface_area
      //
      //
      // percent difference from power diagram
      // ---------------------------------------------
      // prop_map[atom][9] = volume
      // prop_map[atom][10] = overlap_volume
      // prop_map[atom][11] = surface_area
      // prop_map[atom][12] = interfacial_surface_area
      std::map<const Atom *, std::array<double, 13> > prop_map;

      for (Atom_Iter it = model.atoms_begin (); it != model.atoms_end (); ++it)
        {
          std::array<double, 5> awVd_vol{ 0., 0., 0., 0., 0. };
          if (it->residue_name () == "DA" || it->residue_name () == "DT"
              || it->residue_name () == "DG" || it->residue_name () == "DC")
            {
              awVd_vol = find_neighbors (model, *(it), density, options);
            }
          prop_map[&*it][0] = awVd_vol[0];
          prop_map[&*it][1] = awVd_vol[1];
          prop_map[&*it][2] = awVd_vol[2];
          prop_map[&*it][3] = awVd_vol[3];
          prop_map[&*it][4] = awVd_vol[4];
        }

      Rt T;

      power (model, T);

      for (Rt::Vertex_handle vh : T.finite_vertex_handles ())
        {
          std::array<double, 4> power_vol{ 0., 0., 0., 0. };
          if (vh->info ()->residue_name () == "DA"
              || vh->info ()->residue_name () == "DT"
              || vh->info ()->residue_name () == "DG"
              || vh->info ()->residue_name () == "DC")
            {
              power_vol = subdivide (vh, T, options);
            }
          prop_map[vh->info ()][5] = power_vol[0];
          prop_map[vh->info ()][6] = power_vol[1];
          prop_map[vh->info ()][7] = power_vol[2];
          prop_map[vh->info ()][8] = power_vol[3];
        }

      if (options.output_csv)
        {
          csv << "atom,awVd_volume,awVd_overlap_volume,awVd_surface_area,awVd_"
                 "interfacial_surface_area,maximum_gaussian_curvature,power_"
                 "volume,power_overlap_volume,power_surface_area,power_"
                 "interfacial_surface_area,%diff_volume,%diff_overlap_"
                 "volume,%diff_surface_area,%diff_interfacial_surface_area"
              << std::endl;
        }
      for (auto &atom_prop : prop_map)
        {
          auto &prop = atom_prop.second;
          auto atom = *atom_prop.first;
          prop[9] = prop[5] != 0. ? ((prop[0] - prop[5]) / prop[5]) * 100. : 0.;
          prop[10] = prop[6] != 0. ? ((prop[1] - prop[6]) / prop[6]) * 100. : 0.;
          prop[11] = prop[7] != 0. ? ((prop[2] - prop[7]) / prop[7]) * 100. : 0.;
          prop[12] = prop[8] != 0. ? ((prop[3] - prop[8]) / prop[8]) * 100. : 0.;
          if (options.output_csv)
            {
              csv << atom.atom_name () << "_" << atom.residue_name () << "_"
                  << atom.atom_serial_number ();
            }
          for (std::size_t i = 0; i < outs.size (); ++i)
            {
              atom.temperature_factor () = prop[i];
              if (outs[i].is_open ())
                {
                  outs[i] << atom << std::endl;
                }
              if (options.output_csv)
                {
                  csv << "," << prop[i];
                }
            }
          if (options.output_csv)
            {
              csv << std::endl;
            }
        }

      for (auto &elem : outs)
        {
          if (elem.is_open ())
            {
              elem.close ();
            }
        }
      if (csv.is_open ())
        {
          csv.close ();
        }
    }
}
