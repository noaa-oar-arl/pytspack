from numpy.distutils.core import Extension

ext2 = Extension(name='tension',
                 sources=['pytspack/tension.pyf', 'pytspack/tension.f'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name='pytspack',
          version='0.1',
          description="F2PY Users Guide examples",
          author="Barry D. Baker",
          lisense='MIT',
          author_email="barry.baker@noaa.gov",
          source=['pytspack'],
          packages=['pytspack'],
          ext_modules=[ext2]
          )
