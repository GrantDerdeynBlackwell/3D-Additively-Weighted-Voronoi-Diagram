# 3D Additively Weighted Voronoi Diagram
A c++ algorithm to compute the Voronoi diagram of spheres from a molecular structure.

The Voronoi diagram of spheres (sometimes called the additively weighted Voronoi diagram, the Voronoi diagram of balls, or the Apollonius diagram) partitions space around a set of spheres into cells such that each cell contains all the space closer to the surface of a given sphere than any other sphere. It has been postulated that this diagram is better suited to analyze molecular structures than the power diagram. This algorithm is based on the algorithm presented in [[1]](https://doi.org/10.1080/16864360.2016.1273576) and [[2]](https://doi.org/10.1016/j.cag.2019.06.007), though many parts have been adapted or implemented differently.

![awVd_vs_power](https://user-images.githubusercontent.com/104375042/165127476-7d4d1a81-c8b4-427e-8fa5-846fa985b8f3.svg)





## Prereqiuisites

- CGAL >= 5.3
- Eigen >= 3.3
- boost >= 1.71
- ESBTL = 1.0-beta01
- overlap = 0.1.0

## Installation

Collect the dependencies and add them to your CMakeLists.txt. CGAL, Eigen, and boost may be available from your distribution's package manager. An example CmakeLists.txt is included in this repository.
- CGAL: https://www.cgal.org/download.html
- Eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page
- boost: https://www.boost.org/users/download/
- ESBTL: http://esbtl.sourceforge.net/
- overlap: https://github.com/severinstrobl/overlap

Clone this repository. First, generate a makefile with cmake. If there are errors, it is sometimes useful to use cmake-gui.

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

By default, a csv named volumes.csv is output to the current directory. The following flags can be used to change the output behavior:

- ```-vp```: output a pdb where the temperature factor of each atom is the volume of its power cell.
- ```-vv```: output a pdb where the temperature factor of each atom is the volume of its additively weighted cell.
- ```-vd```: output a pdb where the temperature factor of each atom is the percent difference in the volume of its additively weighted cell from the volume of its power cell. (V<sub>aw</sub> - V<sub>p</sub>)/V<sub>p</sub>
- ```-ap```: output a pdb where the temperature factor of each atom is the surface area of its power cell.
- ```-av```: output a pdb where the temperature factor of each atom is the surface area of its additively weighted cell.
- ```-ad```: output a pdb where the temperature factor of each atom is the percent difference in the volume of its additively weighted cell from the volume of its power cell. (A<sub>aw</sub> - A<sub>p</sub>)/A<sub>p</sub>
- ```-k ```: output a pdb where the temperature factor of each atom is the maximum Gaussian curvature of the surface of its additively weighted cell.
- ```-c ```: change the name of the .csv
- ```-mp```: write an OFF mesh for power cell of the molecular structure. Names follow the format power_[ATOM-NAME]_[ATOM-SERIAL#].off
- ```-mv```: write an OFF mesh for additively weighted cell of the molecular structure. Names follow the format awVd_[ATOM-NAME]_[ATOM-SERIAL#].off


The radii associated with each atom can be changed by editing ```data/bondi_classifier.txt```.

## Citations

1.  Zhongyin Hu, Xiang Li, Adarsh Krishnamurthy, Iddo Hanniel & Sara McMains (2017) Voronoi cells of non-general position spheres using the GPU, Computer-Aided Design and Applications, 14:5, 572-581, DOI: 10.1080/16864360.2016.1273576 
2.  Xiang Li, Adarsh Krishnamurthy, Iddo Hanniel, & Sara McMains (2019) Edge topology construction of Voronoi diagrams of spheres in non-general position, Computers & Graphics, 82, 332-342, DOI: 10.1016/j.cag.2019.06.007
