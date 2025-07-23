## Running tests

All tests are stored in the [`test/`](https://github.com/nzfeng/signed-heat-python/tree/main/test) subdirectory, which also contains a few small input files. We use [pytest](https://docs.pytest.org/en/stable/) as a testing framework, which can be installed via `pip install pytest`.

Run the tests with:
```sh
pytest test.py [-q] [-v]
```
where `-v`, `-q` are optional flags that increase/decrease verbosity.
