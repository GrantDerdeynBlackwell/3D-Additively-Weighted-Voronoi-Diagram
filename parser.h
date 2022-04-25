#ifndef G_PARSER
#define G_PARSER

#include <filesystem>

struct Cmd_line_options
{
  bool output_power = false;
  std::filesystem::path power_pdb_fname{ "./power_volumes.pdb" };
  std::filesystem::path power_pdb_overlap_fname{
    "./power_volumes_overlap.pdb"
  };

  bool output_awVd = false;
  std::filesystem::path awVd_pdb_fname{ "./awVd_volumes.pdb" };
  std::filesystem::path awVd_pdb_overlap_fname{ "./awVd_volumes_overlap.pdb" };

  bool output_diff = false;
  std::filesystem::path diff_pdb_fname{ "./diff_volumes.pdb" };
  std::filesystem::path diff_pdb_overlap_fname{ "./diff_volumes_overlap.pdb" };

  bool output_power_areas = false;
  std::filesystem::path pow_area_pdb_fname{ "./power_surface_areas.pdb" };
  std::filesystem::path pow_interface_area_fname{
    "./power_interfacial_areas.pdb"
  };

  bool output_awVd_areas = false;
  std::filesystem::path awVd_area_pdb_fname{ "./awVd_areas.pdb" };
  std::filesystem::path awVd_interface_area_fname{
    "./awVd_interfacial_areas.pdb"
  };

  bool output_areas_diff = false;
  std::filesystem::path diff_area_pdb_fname{ "./diff_areas.pdb" };
  std::filesystem::path diff_area_pdb_interface_fname{
    "./diff_interfacial_areas.pdb"
  };

  bool output_curvature = false;
  std::filesystem::path curvature_pdb_fname{ "./curvatures.pdb" };

  bool mesh_power = false;
  std::filesystem::path power_mesh_dir
      = std::filesystem::exists ("./outputs/") ? "./outputs/" : "./";

  bool mesh_awVd = false;
  std::filesystem::path awVd_mesh_dir
      = std::filesystem::exists ("./outputs/") ? "./outputs/" : "./";

  bool output_csv = true;
  std::filesystem::path csv{ "volumes.csv" };
};

#endif
