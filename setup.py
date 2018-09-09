from setuptools import setup, find_packages

version = open('VERSION').read().strip()

setup(
    name='NavierStokes',
    version=version,
    author='Kyu Mok (Ricky) Kim',
    author_email='rickykim93@hotmail.com',
    url='http://rickykim.net',
    description='Navier Stokes Calculator',
    long_description=open('README.md').read().strip(),
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
    ],
    test_suite='tests',
    entry_points={
        'console_scripts': [
            'packagename = packagename.__main__:main',
        ]
    }
)
