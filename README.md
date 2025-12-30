# pytspack

Python Wrapper around Robert J. Renka's [TSPACK](http://www.netlib.no/netlib/toms/716).

> TSPACK is a curve-fitting package based on exponential tension splines with automatic selection of tension factors.

## Installation

For NOAA Hera users:
```
OPT='-std=c99' pip install https://github.com/noaa-oar-arl/pytspack/archive/master.zip
```
or see the Linux notes [below](#linux).

---

If you are lucky,
```
pip install https://github.com/noaa-oar-arl/pytspack/archive/master.zip
```
will *just work*.

Alternatively, if the above [^b] fails, try
```
pip install https://github.com/noaa-oar-arl/pytspack/archive/master.zip --no-use-pep517 --no-deps
```
(you must already have `numpy` installed and currently must have `setuptools` earlier than v65).

Otherwise, you can clone the repo and try to build the extension module
using [`f2py`](https://numpy.org/doc/stable/f2py/index.html) manually...


[^b]: [Build-time dependency specification](https://pip.pypa.io/en/stable/reference/build-system/pyproject-toml/#build-time-dependencies) via `pyproject.toml` file.

### Windows

On Windows, with a [GNU Fortran via MSYS2](https://numpy.org/doc/stable/f2py/windows/msys2.html)
setup, try the following:

1. Clone the repo
   ```powershell
   git clone https://github.com/noaa-oar-arl/pytspack
   cd pytspack
   ```

2. Build the extension, using one of these options

   Using the `setup.py`:
   ```
   python setup.py build_ext --fcompiler=gnu95 --compiler=mingw32 --build-lib=./pytspack
   ```

   OR

   Using `f2py` directly:
   ```powershell
   cd pytspack
   python -m numpy.f2py -c --fcompiler=gnu95 --compiler=mingw32 -m tspack tspack.f
   cd ..
   ```

3. Link pytspack to your active Python environment.
   ```powershell
   pip install -e . --no-use-pep517
   ```

### Linux

We have seen some issues where some `gcc`s do not want to compile the `fortranobject.c` [^a]

Follow the same general steps as above, but don't use `--compiler=mingw32`.

On NOAA Hera, the default GCC is currently v4.
Use `module load gnu` to get newer before attempting to install pytspack,
or use
```bash
OPT='-std=c99' python -m numpy.f2py -c -m tspack tspack.f
```
to configure `gcc`.

`OPT='-std=c99'` can also be prepended to [the `pip install` above](#installation),
enabling an installation without cloning the repo.


[^a]: As mentioned [here](https://mfix.netl.doe.gov/forum/t/strange-build-error-in-mfix-21-4/3923/3), for example.

## More information

For more information on TSPACK, see [the open-access paper](https://dl.acm.org/doi/10.1145/151271.151277) (Renka, 1993) or [the netlib page](https://web.archive.org/web/20220701070001/http://netlib.no/netlib/toms/716).

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
