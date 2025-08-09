# signedheat3d (Python package)

## Documentation lives at [https://nzfeng.github.io/signed-heat-3d/py](https://nzfeng.github.io/signed-heat-3d/py)

![teaser image](https://github.com/nzfeng/signed-heat-3d/blob/main/media/teaser.png)

A Python library implementing the [Signed Heat Method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) for computing robust signed distance fields (SDFs) to polygon meshes and point clouds in 3D.

* The original C++ code lives at [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d).
* If you're interested in using the Signed Heat Method *on* 2D surface domains, rather than in 3D Euclidean space, check out the signed heat method implementation in [`potpourri3d`](https://github.com/nmwsharp/potpourri3d). (The overall organization of this repository was inspired by that of `potpourri3d`!)

## Installation

The recommended way to install `signedheat3d` is via PyPI:

```
pip install signedheat3d
```
You can also clone the repository and install it from source:
```
git clone --recurse-submodules https://github.com/nzfeng/signed-heat-python.git
cd signed-heat-python
pip install .
```
If you do not clone recursively, some submodules or sub-submodules will not clone. Initialize/update these submodules by running `git submodule update --init --recursive` or `git submodule update --recursive`.

<!-- If you are a developer, then it is usually much more efficient to install the build dependencies in your environment once and use the following command that avoids creating a new virtual environment at every compilation:

```
pip install --no-build-isolation -ve .
``` -->

### Dependencies

This project has the following submodules, which should have been installed in the previous step.
* [nanobind](https://nanobind.readthedocs.io/en/latest/)
* [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d)

You may need to install Boost since `signed-heat-3d` depends on [`amgcl`](https://github.com/ddemidov/amgcl), which depends on Boost. Boost can be installed on macOS using `brew install boost`, and the necessary modules on Ubuntu using
```
sudo apt-get -y update
sudo apt-get -y install libboost-dev libboost-test-dev libboost-program-options-dev libboost-serialization-dev
```
Windows users should follow the instructions on the [Boost website](https://www.boost.org/releases/latest/).

## Demo program

This repository also contains a demo Python program at `test/demo.py`, using [Polyscope](https://github.com/nmwsharp/polyscope-py) for visualization. The demo program at `test/demo.py` uses the following Python libraries, which can each be installed via `pip install`:
* [NumPy](https://numpy.org/)
* [polyscope](https://polyscope.run/py/)
<!-- * [mypy](https://www.mypy-lang.org/) (assuming Python 3.8+) -->
To run the demo program, pip-install the package using the instructions described [above](#installation). Then `cd` into the top level of the directory, and run
```

python3 test/demo.py path/to/mesh/or/pointcloud
```

### Input / Output

Input / output meshes can be any one of the following [types](https://geometry-central.net/surface/utilities/io/) ("Supported file types"), including OBJ, PLY, STL, and OFF.

Point clouds are currently assumed to have file extension `.pc` and consist of newline-separated 3D point positions (denoted by leading char `v`) and point normal vectors (denoted by leading char `vn`).

### Command line arguments

In addition to the mesh file, you can pass several flags.

|flag | purpose|
| ------------- |-------------|
|`--g`, `--grid`| Solve on a background grid. By default, the domain will be discretized as a tet mesh. |
|`--v`, `--verbose`| Verbose output. Off by default.|
|`--s`| Controls the tet/grid spacing proportional to $2^{-h}$, with larger values indicating more refinement. Default value is 0.|
|`--l`, `--headless`| Don't use the GUI, and automatically solve for & export the generalized SDF.|

<!-- ## TODOs

* Contouring slower than in [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d), because data is being passed by value with each call to the Python-bound functions
* More precise level set constraints for grid solves
* Isoline rendering for volume meshes is [not yet bound in Polyscope](https://github.com/nmwsharp/polyscope-py/issues/36); for now, SDFs can be rendered with isobands via the GUI only.
* Handle more input file formats, via extra Python bindings to [geometry-central](https://geometry-central.net/)'s IO functions.
 -->