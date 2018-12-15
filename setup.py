#!/usr/bin/python

from setuptools import setup
from setuptools import find_packages
from os.path import join as opj
from os.path import dirname

def get_version():
    """Load version only
    """
    with open(opj(dirname(__file__), 'multimatch', '__init__.py')) as f:
        version_lines = list(filter(lambda x: x.startswith('__version__'), f))
    assert (len(version_lines) == 1)
    return version_lines[0].split('=')[1].strip(" '\"\t\n")

# extension version
version = get_version()



README = opj(dirname(__file__), 'README.md')
try:
    import pypandoc
    long_description = pypandoc.convert(README, 'rst')
except (ImportError, OSError) as exc:
    print(
        "WARNING: pypandoc failed to import or threw an error while converting"
        " README.md to RST: %r  .md version will be used as is" %exc
    )
    long_description = open(README).read

# Metadata
setup(
    name='multimatch',
    version=version,
    description='Multidimensional scanpath comparison',
    long_description=long_description,
    author='Adina Wagner',
    author_email='adina.wagner@t-online.de',
    url='https://github.com/AdinaWagner/multimatch',
    packages=[pkg for pkg in find_packages('.') if pkg.startswith('multimatch')],
    install_requires=[
        'numpy',
        'pandas',
        #'math',
        #'os',
        #'bisect,       #all of these requirements fail with testpypi?
    ],
    extras_require={
        'devel-docs': [
            # for converting README.md -> .rst for long description
            'pypandoc',
        ]},
    entry_points={
        'console_scripts': [
            'multimatch=multimatch.multimatch:main',
            'multimatch_forrest=multimatch.multimatch_forrest:main'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

