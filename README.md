# 3D Additively Weighted Voronoi Diagram
A c++ algorithm to compute the Voronoi diagram of spheres from a molecular structure.

The Voronoi diagram of spheres (sometimes called the additively weighted Voronoi diagram, the Voronoi diagram of balls, or the Apollonius diagram) partitions space around a set of spheres into cells such that each cell contains all the space closer to the surface of a given sphere than any other sphere. It has been postulated that this diagram is better suited to analyze molecular structures than the power diagram. This algorithm is based on the algorithm presented in [[1]](https://doi.org/10.1080/16864360.2016.1273576) and [[2]](https://doi.org/10.1016/j.cag.2019.06.007), though many parts have been adapted or implemented differently.

![A diagram illustrating the difference between a power cell and an additively weighted cell](./images/power_vs_voronoi.svg?raw=true "Power cell vs additively weighted cell")


## Prereqiuisites

- CGAL >= 5.3
- Eigen >= 3.3
- boost >= 1.71
- ESBTL = 1.0-beta01
- overlap = 0.1.0

## Installation

### Distro specific installation
Copy and paste the following into a terminal to easily collect each dependency and compile.

<details>
  <summary>Ubuntu</summary>
  
  ```bash
 sudo apt install libeigen3-dev libboost-all-dev
 git clone https://github.com/severinstrobl/overlap.git
 wget https://github.com/CGAL/cgal/archive/refs/tags/v5.4.tar.gz
 tar -xvf v5.4.tar.gz
 wget https://versaweb.dl.sourceforge.net/project/esbtl/ESBTL-1.0-beta01.tar.bz2
 tar -xvf ESBTL-1.0-beta01.tar.bz2
 git clone https://github.com/GrantDerdeynBlackwell/3D-Additively-Weighted-Voronoi-Diagram.git
 cd 3D-Additively-Weighted-Voronoi-Diagram
 cmake -DCGAL_DIR="../cgal-5.4" -DCMAKE_BUILD_TYPE="Release" .
 make
 ```
</details>

### Manual Installation

Collect the dependencies and add them to your CMakeLists.txt. CGAL, Eigen, and boost may be available from your distribution's package manager. An example CmakeLists.txt is included in this repository.
- CGAL: https://www.cgal.org/download.html
- Eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page
- boost: https://www.boost.org/users/download/
- ESBTL: http://esbtl.sourceforge.net/
- overlap: https://github.com/severinstrobl/overlap

Clone this repository. First, generate a makefile with cmake. If there are errors, it is sometimes useful to use cmake-gui. Most of the dependencies can be automatically found with `cmake` or at worst by passing the appropriate include directory to `cmake-gui`. However, you must explicitly add the ESBTL include directory to `CMakeLists.txt`.

Once you have collected all the dependencies and added them to CMake:

```bash
cmake .
```
or
```bash
cmake-gui .
```


Compile with a C++14 compatible compiler such as GCC or Clang.

```bash
make
```
    
## Usage

```bash
./vdos FILENAME SUBDIVISIONS [OPTIONS]
```

FILENAME is the filename of the pdb you wish to analyze.

SUBDIVISIONS controls the density of the ray sampling. A higher number will produce a better approximation. We have found that 2 subdivions is a good trade-off between speed and accuracy.

### Include/Exclude

By default, vdos only computes the power or awVd cells of atoms whose belong to the following residues:

    DA,  DG,  DT,  DC,  ALA, ARG, ASN, ASP, ASX,
    
    CYS, CYX, GLU, GLN, GLX, GLY, HIS, ILE, LEU,
    
    LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
    
To change this behavior, use the ```--ra``` flag to add residue names or the ```--rs``` to subtract residue names. For example, to include residues name NA, CL, K, & BR, and exlcude residues named DC & DT:

```bash
./vdos example.pdb 2 --ra NA CL K BR --rs DC DT
```

### Radii

By default, vdos uses the file ```data/bondi_classifier.txt``` to associate radii to each atom. To change the radii, you can either edit ```data/bondi_classifier``` or you can change the b-factor in your input pdb to the radii. If you use the latter, you must pass the flag ```--radii_from_b```.

### Output

By default, a csv named volumes.csv is output to the current directory. The following flags can be used to change the output behavior:

- ```--vp```: output a pdb where the temperature factor of each atom is the volume of its power cell. Default: power_volumes.pdb
- ```--vv```: output a pdb where the temperature factor of each atom is the volume of its additively weighted cell. Default: awVd_volumes.pdb
- ```--vd```: output a pdb where the temperature factor of each atom is the percent difference in the volume of its additively weighted cell from the volume of its power cell. (V<sub>aw</sub> - V<sub>p</sub>)/V<sub>p</sub> Default: diff_volumes.pdb
- ```--op```: output a pdb where the temperature factor of each atom is the volume of its power cell occupied by an atom. Default: power_overlap_volumes.pdb
- ```--ov```: output a pdb where the temperature factor of each atom is the volume of its awVd cell occupied by an atom. Default: awVd_overlap_volumes.pdb
- ```--od```: output a pdb where the temperature factor of each atom is the percent difference in the volume of its power cell occupied by an atom from the volume of its awVd cell occupied by an atom. Default: diff_overlap_volumes.pdb
- ```--ap```: output a pdb where the temperature factor of each atom is the surface area of its power cell. Default: power_surface_areas.pdb
- ```--av```: output a pdb where the temperature factor of each atom is the surface area of its additively weighted cell. Default: awVd_surface_areas.pdb
- ```--ad```: output a pdb where the temperature factor of each atom is the percent difference in the volume of its additively weighted cell from the volume of its power cell. (A<sub>aw</sub> - A<sub>p</sub>)/A<sub>p</sub> Default: diff_surface_areas.pdb
- ```--ip```: output a pdb where the temperature factor of each atom is the sum of the surface area of each face of the power cell at an interface. That is, each face between an atom in the residue_list and an atom not in the residue_list. Default: power_interface_areas.pdb
- ```--iv```: output a pdb where the temperature factor of each atom is the sum of the surface area of each face of the awVd cell at an interface. That is, each face between an atom in the residue_list and an atom not in the residue_list. Default: awVd_interface_areas.pdb
- ```--id```: output a pdb where the temperature factor of each atom is percent difference of the interfacial area of its awVd cell from its power cell. Default: diff_interface_areas.pdb
- ```-k ```: output a pdb where the temperature factor of each atom is the maximum Gaussian curvature of the surface of its additively weighted cell. Default: max_curvature.pdb
- ```-c ```: change the name of the .csv
- ```--verts```: output the vertices of the awVd to stdout in the format {a1, a2, a3, a4, x, y, z}
- ```--mp```: write an OFF mesh for power cell of the molecular structure. Names follow the format power_[ATOM-NAME]_[ATOM-SERIAL#].off
- ```--mv```: write an OFF mesh for additively weighted cell of the molecular structure. Names follow the format awVd_[ATOM-NAME]_[ATOM-SERIAL#].off

## Examples

| ![Some representative cells of atoms. Faces are colored randomly for clarity.](./images/cell_picutres.jpeg?raw=true "Representative additively weighted cells") |
|:--:| 
| *Additively weighted Voronoi cells of individul atoms. Top to bottom: hydrogen, carbon, nitrogen, oxygen, phosphorus. Left to right: Atoms contributing to the cell, atoms and the cell, atoms with radii reduced by a constant factor, additively weighted Voronoi cell.* |

| ![The additively weighted Voronoi shell of 1BNA. Faces are colored by the element which the face contacts.](./images/DNA_example.jpeg?raw=true "The additively weighted Voronoi shell of 1BNA") |
|:--:| 
| *Additively weighted Voronoi cells of 1BNA atoms. Left: Cartoon of DNA structure. Waters are represented as bent lines. Middle-top: Atomic sphere representation of DNA. Middle-bottom: Additively weighted Voronoi diagram of DNA. Right-top: Atomic sphere representation of DNA. Waters are represented as translucent spheres. Right-bottom: Additively weighted Voronoi diagram of DNA with waters as translucent spheres.* |

## Citations

1.  Zhongyin Hu, Xiang Li, Adarsh Krishnamurthy, Iddo Hanniel & Sara McMains (2017) Voronoi cells of non-general position spheres using the GPU, Computer-Aided Design and Applications, 14:5, 572-581, DOI: 10.1080/16864360.2016.1273576 
2.  Xiang Li, Adarsh Krishnamurthy, Iddo Hanniel, & Sara McMains (2019) Edge topology construction of Voronoi diagrams of spheres in non-general position, Computers & Graphics, 82, 332-342, DOI: 10.1016/j.cag.2019.06.007
