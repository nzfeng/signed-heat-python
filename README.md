# signed-heat-python

A Python library implementing the [Signed Heat Method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) for computing robust signed distance fields (SDFs) to polygon meshes and point clouds in 3D.

* [signed-heat-demo](https://github.com/nzfeng/signed-heat-demo)
* [signed-heat-3d](https://github.com/nzfeng/signed-heat-3d)
* For more geometry processing tools on surface meshes in Python, check out [`potpourri3d`](https://github.com/nmwsharp/potpourri3d). (The overall organization of this repository was inspired by that of `potpourri3d`!)

## Installation

```
git clone --recursive https://github.com/nzfeng/signed-heat-python.git
```

If you do not clone recursively, some submodules or sub-submodules will not clone. Initialize/update these submodules by running `git submodule update --init --recursive` or `git submodule update --recursive`.

## Getting started

To run the demo program, 

```
cd signed-heat-python
mkdir build && cd build
cmake .. && make -j
python3 ../demo.py
```

## Documentation

### Input / Output

For simplicity, input / output mesh files are assumed to be OBJ format, though this could be amended with an extra Python binding to [geometry-central](https://geometry-central.net/)'s IO functions.

### Signed heat method
