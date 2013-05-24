from distutils.core import setup, Extension
# Run with python setup.py build_ext --inplace

setup(
    name='Needleman-Wunsch', version='1.0',
    ext_modules=[Extension('nwmodule', ['nwmodule.c'],
                           extra_compile_args=['-std=c99'])]
    )
