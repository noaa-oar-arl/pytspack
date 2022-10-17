from numpy.distutils.core import Extension

ext = Extension(
    name="tspack",
    sources=["pytspack/tspack.pyf", "pytspack/tspack.f"],
)

if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(
        name="pytspack",
        version="0.1.2",
        description="Wrapper around Robert J. Renka's fortran TSPACK: Tension Spline Curve Fitting Package",
        author="Barry D. Baker",
        license="MIT",
        author_email="barry.baker@noaa.gov",
        packages=["pytspack"],
        ext_modules=[ext],
        install_requires=["numpy"],
    )
