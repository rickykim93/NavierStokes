# NavierStokes
Pre v.1.0.0 Release  - Still in development
This is a Navier Stokes calculator using FVM 
(Finite Volume Method) which is used
in computational fluid dynamics.
This is the same formula that ANSYS Fluent uses.

## Getting Started
### Prerequisites
* `Numpy`
* `Scipy`
* `MatplotLib`

### Installing
```commandline
pip install NavierStokes
```
To upgrade:
```commandline
pip install NavierStokes -U
```

## Running the tests
Requires `pytest`
```commandline
python -m unittest tests.test_TDMAsolver
```

## Documentation
The documentation for:
* [Command Line Functions](https://github.com/rickykim93/NavierStokes/blob/master/NavierStokes/README.md)
* [TDMA](https://rickykim.net/NSDoc#tdma) - LinkWIP
* [Diffusion](https://rickykim.net/NSDoc#diffusion) - LinkWIP
* [Advection](https://rickykim.net/NSDoc#advection) - LinkWIP
* [SIMPLE](https://rickykim.net/NSDoc#simple) - LinkWIP

## Example - See it in Action!
[This page](https://rickykim.net/NavierStokes) uses this module. - LinkWIP

## Contributing
All contributions, bug reports, bug fixes, documentation improvements, enhancements and ideas are welcome.

Report all bugs and issues to [this page](https://github.com/rickykim93/NavierStokes/issues)

## Versioning
I use [SemVer](https://semver.org) for versioning.
For the versions available, see the [tags on this repository](https://github.com/rickykim93/NavierStokes/releases)

## Authors
* Kyu Mok (Ricky) Kim - Initial work - [website](https://rickykim.net), [github](https://github.com/rickykim93)

## License
This project is licensed under the MIT License - see the [License](https://github.com/rickykim93/NavierStokes/blob/master/LICENSE)
file for details.

## Acknowledgements
* [MIE 1210: Computational Fluid Mechanics](https://rickykim.net) - LinkWIP Course from University of Toronto
for teaching me Navier Stokes equations and how to implement them
