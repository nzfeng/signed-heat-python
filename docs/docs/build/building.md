# Building

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

#### Example

For a simple example project using [Polyscope](https://polyscope.run) for visualization, see the sample Python project [here](https://github.com/nzfeng/signed-heat-python/blob/main/test/demo.py).
