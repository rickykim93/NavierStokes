from setuptools import setup, find_packages
from NavierStokes import __version__
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='NavierStokes',
    version=__version__,
    author='Kyu Mok (Ricky) Kim',
    author_email='rickykim93@hotmail.com',
    url='https://github.com/rickykim93/NavierStokes',
    description='Navier Stokes Calculator',
    long_description=open(path.join(here, 'README.md'), encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
    ],
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'NavierStokes = NavierStokes.__main__:main',
        ]
    }
)
