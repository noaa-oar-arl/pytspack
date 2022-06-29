# pytspack

Python Wrapper around Robert J. Renka's [TSPACK](http://www.netlib.no/netlib/toms/716).
TSPACK

> TSPACK is a curve-fitting package based on exponential tension splines with automatic selection of tension factors.

## Installation

If you are lucky,
```
pip install https://github.com/noaa-oar-arl/pytspack/archive/master.zip --no-deps
```
will work.

On Windows, a [GNU Fortran via MSYS2](https://numpy.org/doc/stable/f2py/windows/msys2.html)
setup is expected.

Otherwise, you can clone the repo and try to build the extension module
using [F2PY](https://numpy.org/doc/stable/f2py/index.html) manually.

## Citation

For more information on TSPACK, see [the open-access paper](https://dl.acm.org/doi/10.1145/151271.151277) (Renka, 1993) or [the netlib page](http://www.netlib.no/netlib/toms/716).
