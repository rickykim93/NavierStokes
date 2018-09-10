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
    url='http://rickykim.net',
    description='Navier Stokes Calculator',
    long_description=open(path.join(here, 'README.md'), encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
    ],
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'NavierStokes = NavierStokes.__main__:main',
        ]
    }
)
