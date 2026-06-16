# Python Wrapper

The Python wrapper exposes the Splines C++ library through `pybind11`.

These instructions are for the `src_python_interface` directory and use
**Conda** for Python, CMake, compiler toolchain, and `pybind11`.

## Requirements

Install Miniconda or Miniforge. On macOS, for example:

```bash
brew install --cask miniconda
conda init "$(basename "${SHELL}")"
```

Then restart the terminal.

## Create the Conda environment

From `src_python_interface`:

```bash
conda env create -f environment.yml
conda activate splines-py
```

Check that the environment is coherent:

```bash
python --version
cmake --version
python -c "import pybind11; print(pybind11.get_cmake_dir())"
```

## Configure and build

The wrapper CMakeLists uses modern `FindPython`, so the important variable is
`Python_EXECUTABLE`.

From `src_python_interface`:

```bash
rm -rf build

cmake -S . -B build \
  -DPython_EXECUTABLE:FILEPATH="${CONDA_PREFIX}/bin/python"

cmake --build build --parallel
cmake --install build
```

The installed module is written to:

```text
../lib/lib/Splines.{ext}
```

On macOS, for example:

```text
../lib/lib/Splines.cpython-312-darwin.so
```

## If CMake picks the wrong Python or cannot find pybind11

On some systems CMake may resolve the base Conda interpreter instead of the
active environment. In that case configure explicitly with:

```bash
cmake -S . -B build \
  -DPython_EXECUTABLE:FILEPATH="${CONDA_PREFIX}/bin/python" \
  -DPython_ROOT_DIR:PATH="${CONDA_PREFIX}" \
  -DPython_FIND_VIRTUALENV=ONLY \
  -Dpybind11_DIR:PATH="$(python -c 'import pybind11; print(pybind11.get_cmake_dir())')"
```

## Test the wrapper

After the build:

```bash
ctest --test-dir build --output-on-failure
```

Or run the smoke test directly against the installed module:

```bash
python test_splines.py --module-dir ./../lib/lib
```

You can also test the module directly from the build directory before install:

```bash
python test_splines.py --module-dir ./build
```

## Usage

In Python:

```python
import sys
sys.path.insert(0, "../lib/lib")

import Splines

help(Splines)
```

The wrapper follows the current C++ API. For example, the monotone cubic
variant exposed in Python is `VanLeerSpline`.

## Notes

- `libSplines_*` depends on `libquarticRootsFlocke_*`, and the wrapper links it automatically through CMake.
- If you switch Python environments, remove `build` and reconfigure.
- Do not install `pybind11` manually from GitHub unless you explicitly need a custom upstream revision.
