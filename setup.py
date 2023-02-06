#!/usr/bin/env python3
import setuptools
from setuptools.extension import Extension

try:
    # Build-time dependencies... specified in pyproject.toml but being careful anyway.
    import Cython.Build
    import numpy
except ImportError:
    print("Please install 'Cython' and 'numpy' before installing this package!")
    raise


with open("README.md", "rt", encoding="UTF-8") as fh:
    long_description = fh.read()

version = {}
with open("ezstructure/__init__.py", "rt", encoding="UTF-8") as fh:
    exec(fh.read(), version)

numpy_include = numpy.get_include()


setuptools.setup(
    name="ezstructure",
    version=version["__version__"],
    description="EasyStructure",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Bug Tracker": "https://github.com/Jerrythafast/ezstructure/issues",
        "Source Code": "https://github.com/Jerrythafast/ezstructure",
    },
    author="Jerry Hoogenboom",
    author_email="jerryhoogenboom@outlook.com",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords="bioinformatics forensics structure ancestry",
    packages=setuptools.find_packages(),
    ext_modules=[
        Extension(
            "ezstructure.vars.utils",
            sources=["ezstructure/vars/utils.pyx"],
            include_dirs=[numpy_include]),
        Extension(
            "ezstructure.vars.admixprop",
            sources=["ezstructure/vars/admixprop.pyx", "ezstructure/vars/C_admixprop.c"],
            include_dirs=[numpy_include]),
        Extension(
            "ezstructure.vars.allelefreq",
            sources=["ezstructure/vars/allelefreq.pyx", "ezstructure/vars/C_allelefreq.c"],
            #libraries=["gsl"],
            extra_compile_args=["-O3"],
            include_dirs=[numpy_include]),
        Extension(
            "ezstructure.vars.marglikehood",
            sources=["ezstructure/vars/marglikehood.pyx", "ezstructure/vars/C_marglikehood.c"],
            include_dirs=[numpy_include]),
        Extension(
            "ezstructure.fastStructure",
            sources=["ezstructure/fastStructure.pyx"],
            include_dirs=[numpy_include, 'ezstructure/vars/'])],
    cmdclass={"build_ext": Cython.Build.build_ext},
    package_data={
        "ezstructure": ["vars/*.pxd", "vars/*.h"]
    },
    install_requires=["numpy", "scipy"],
    python_requires=">=3.5",
    entry_points={
        "console_scripts": [
            "ezstructure=ezstructure.structure:main",
            "ezdistruct=ezstructure.distruct:main",
            "ezstructure_choosek=ezstructure.chooseK:main"
        ]
    },
    zip_safe=False  # Needed for Cython, according to Cython docs.
)
