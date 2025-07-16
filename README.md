# signed_heat_method

![teaser image](https://github.com/nzfeng/signed-heat-3d/blob/main/media/teaser.png)

A Python library implementing the [Signed Heat Method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) for computing robust signed distance fields (SDFs) to polygon meshes and point clouds in 3D.

* The original C++ code lives at [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d).
* If you're interested in using the Signed Heat Method *on* 2D surface domains, rather than in 3D Euclidean space, check out [signed-heat-demo](https://github.com/nzfeng/signed-heat-demo).
* For more geometry processing tools on surface meshes in Python, check out [`potpourri3d`](https://github.com/nmwsharp/potpourri3d). (The overall organization of this repository was inspired by that of `potpourri3d`!)

## Installation

The recommended way to install `signed_heat_method` is via PyPI:

```
pip install signed_heat_method
```
You can also clone the repository and install it from source:
```
git clone --recursive https://github.com/nzfeng/signed-heat-python.git
cd signed-heat-python
git submodule update --init --recursive
pip install .
```
If you do not clone recursively, some submodules or sub-submodules will not clone. Initialize/update these submodules by running `git submodule update --init --recursive` or `git submodule update --recursive`.

If you are a developer, then it is usually much more efficient to install the build dependencies in your environment once and use the following command that avoids creating a new virtual environment at every compilation:

```
pip install --no-build-isolation -ve .
```

### Dependencies

This project has the following submodules, which should have been installed in the previous step.
* [nanobind](https://nanobind.readthedocs.io/en/latest/)
* [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d)

You may need to install Boost since `signed-heat-3d` depends on [`amgcl`](https://github.com/ddemidov/amgcl), which depends on Boost (`brew install boost` on macOS).

The demo program at `test/demo.py` uses the following Python libraries, which can each be installed via `pip install`:
* [NumPy](https://numpy.org/)
* [polyscope](https://polyscope.run/py/)
* [mypy](https://www.mypy-lang.org/) (assuming Python 3.8+)

## Getting started

Documentation is [below](#documentation). This repository also contains a demo Python program at `test/demo.py`, using [Polyscope](https://github.com/nmwsharp/polyscope-py) for visualization. To run the demo program, pip-install the package using the instructions described [above](#installation).

Then `cd` into the top level of the directory, and run
```

python3 test/demo.py path/to/mesh/or/pointcloud
```

## Documentation

### Input / Output

For simplicity, input / output mesh files are assumed to be OBJ format, though this could be amended with an extra Python binding to [geometry-central](https://geometry-central.net/)'s IO functions.

Point clouds are currently assumed to have file extension `.pc` and consist of newline-separated 3D point positions (denoted by leading char `v`) and point normal vectors (denoted by leading char `vn`).

### Command line arguments

In addition to the mesh file, you can pass several flags.

|flag | purpose|
| ------------- |-------------|
|`--g`, `--grid`| Solve on a background grid. By default, the domain will be discretized as a tet mesh. |
|`--v`, `--verbose`| Verbose output. Off by default.|
|`--s`| Controls the tet/grid spacing proportional to $2^{-h}$, with larger values indicating more refinement. Default value is 0.|
|`--l`, `--headless`| Don't use the GUI, and automatically solve for & export the generalized SDF.|

To improve performance, operators and spatial discretizations are only built as necessary, and re-used in future computations if the underlying discretization hasn't changed. This means future computations can be significantly faster than the initial solve (which includes, for example, tet mesh construction and matrix factorization.)

### Signed distance to a mesh

```python
import signed_heat_method as shm

V, F = # your mesh

# Initalize tet mesh solver
tet_solver = shm.SignedHeatTetSolver(verbose=True) # print all messages

# Solve!
sdf_tets = tet_solver.compute_distance_to_mesh(V, F)

# Solve on a regular grid instead.
grid_solver = shm.SignedHeatGridSolver(verbose=True)
sdf_grid = grid_solver.compute_distance_to_mesh(V, F)

```

- `SignedHeatTetSolver.compute_distance_to_mesh(V, F, options)`
- `SignedHeatGridSolver.compute_distance_to_mesh(V, F, options)`
  - `V` a |V| x 3 NumPy array of 3D vertex positions, where |V| = number of vertices in the source mesh.
  - `F` a list of lists; each sublist represents a polygonal mesh face of arbitrary degree (using 0-indexed vertices).

`options` is a dictionary with the following possible entries:
  - `level_set_constraint` Whether to apply level set constraints, with options "ZeroSet", "None", "Multiple". Generally set to "ZeroSet" (set by default).
  - `t_coef` Sets the time used for short-time heat flow. Default value is 1. Generally don't change this.
  - `bbox_min` The 3D position of the minimum corner of the computational domain, which is assumed to be an axis-aligned rectangular prism. If not specified, the size of the domain will be automatically computed so as to encompass the input source geometry.
  - `bbox_max` The 3D position of the maximum corner of the computational domain.
  - `resolution` 3D vector specifying the tet/grid spacing, with larger values indicating more refinement. If solving on a grid, this corresponds to the number of nodes along each dimension. Default values are $2^{5}$.
  - `rebuild` If `True`, (re)build the underlying tet mesh or grid domain.

### Signed distance to a point cloud

```python
import signed_heat_method as shm

P, N = # your point cloud, with normals

# Initalize tet mesh solver
tet_solver = shm.SignedHeatTetSolver(verbose=True) # print all messages

# Solve!
sdf_tets = tet_solver.compute_distance_to_point_cloud(P, N)

# Solve on a regular grid instead.
grid_solver = shm.SignedHeatGridSolver(verbose=True)
sdf_grid = grid_solver.compute_distance_to_point_cloud(P, N)

```

- `SignedHeatTetSolver.compute_distance_to_point_cloud(P, N, options)`
- `SignedHeatGridSolver.compute_distance_to_point_cloud(P, N, options)`
  - `P` a |P| x 3 NumPy array of 3D point positions, where |P| = number of points in the source point cloud.
  - `N` a |P| x 3 NumPy array of 3D normal vectors, where |P| = number of points in the source point cloud.

`options` is a dictionary with the following possible entries:
  - `t_coef` Sets the time used for short-time heat flow. Default value is 1. Generally don't change this.
  - `bbox_min` The 3D position of the minimum corner of the computational domain, which is assumed to be an axis-aligned rectangular prism. If not specified, the size of the domain will be automatically computed so as to encompass the input source geometry.
  - `bbox_max` The 3D position of the maximum corner of the computational domain.
  - `resolution` 3D vector specifying the tet/grid spacing, with larger values indicating more refinement. If solving on a grid, this corresponds to the number of nodes along each dimension. Default values are $2^{5}$.
  - `rebuild` If `True`, (re)build the underlying tet mesh or grid domain.

### Auxiliary functions

- `SignedHeatTetSolver.get_vertices()` returns an nVertices x 3 array representing the vertex locations of the underlying tet mesh domain.
- `SignedHeatTetSolver.get_tets()` returns an nTets x 4 array representing the tetrahedra of the underlying tet mesh domain, where each tetrahedra is given by four vertex indices (0-indexed).
- `SignedHeatTetSolver.isosurface(f, isovalue)` contours a scalar function defined on the tet mesh, given by the input vector `f`, according to the isovalue `isovalue`. Returns the tuple `(vertices, faces)` defining the polygon mesh of the resulting isosurface.

- `SignedHeatGridSolver.get_grid_resolution()` returns a length-3 array giving the number of cells of the background grid along the x-, y-, and z-axes, respectively.
- `SignedHeatGridSolver.get_bbox()` returns the tuple `(bbox_min, bbox_min)`, where `bbox_min` and `bbox_min` are the 3D positions of the minimal/maximal node corners of the grid.
- `SignedHeatGridSolver.to_grid_array(f, isovalue)`  converts the input vector `f` representing a scalar function on the background grid, into a NumPy array of shape `(dim_x, dim_y, dim_z)`, where `dim_[xyz]` gives the number of grid nodes along each axis.

## TODOs

* Testing, CI
* Python package release
* Contouring much slower than in [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d), because data is being passed by value with each call to the Python-bound functions
* More precise level set constraints for grid solves
* Isoline rendering for volume meshes is [not yet bound in Polyscope](https://github.com/nmwsharp/polyscope-py/issues/36); for now, SDFs can be rendered with isobands via the GUI only.
* Handle more input file formats, via extra Python bindings to [geometry-central](https://geometry-central.net/)'s IO functions.
