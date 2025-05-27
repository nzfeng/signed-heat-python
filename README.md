# signed-heat-python

A Python library implementing the [Signed Heat Method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) for computing robust signed distance fields (SDFs) to polygon meshes and point clouds in 3D.

* [](https://github.com/nzfeng/signed-heat-demo)
* [](https://github.com/nzfeng/signed-heat-3d)
* [](https://geometry-central.net/)
* For more geometry processing tools on surface meshes in Python, check out [`potpourri3d`](http://geometry-central.net/). (The overall organization of this repository was inspired by that of `potpourri3d`!)

## Installation

## Getting started

To run the demo program, 

```
mkdir build && cd build
cmake .. && make -j
python3 ../demo.py
```

## Documentation

- [Input / Output](#input--output)
- [Mesh basic utilities](#mesh-basic-utilities)

### Input / Output

Read/write meshes and point clouds from some common formats.
