#!/usr/bin/python

from setuptools import setup
from setuptools import find_packages
from os.path import join as opj
from os.path import dirname

def get_version():
    """Load version only
    """
    with open(opj(dirname(__file__), 'MultiMatch_pure', '__init__.py')) as f:
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
    url='https://github.com/AdinaWagner/MultiMatch',
    packages=[pkg for pkg in find_packages('.') if pkg.startswith('MultiMatch')],
    install_requires=[
        'numpy',
        'pandas',
        'math',
        'bisect',
        'os',
    ],
    extras_require={
        'devel-docs': [
            # for converting README.md -> .rst for long description
            'pypandoc',
        ]},
    entry_points={
        'console_scripts': [
            'multimatch=MultiMatch.multimatch:main',
            'multimatch_forrest=MultiMatch.multimatch_forrest:main'],
    }
)

