# Grd2Shp_Xagg

[![PyPI](https://img.shields.io/pypi/v/grd2shp_xagg.svg)](https://pypi.org/project/grd2shp_xagg/)
[![Status](https://img.shields.io/pypi/status/grd2shp_xagg.svg)](https://pypi.org/project/grd2shp_xagg/)
[![Python
Version](https://img.shields.io/pypi/pyversions/grd2shp_xagg)](https://pypi.org/project/grd2shp_xagg)
[![License](https://img.shields.io/pypi/l/grd2shp_xagg)](https://opensource.org/licenses/MIT)

[![Read the documentation at
https://grd2shp_xagg.readthedocs.io/](https://img.shields.io/readthedocs/grd2shp_xagg/latest.svg?label=Read%20the%20Docs)](https://grd2shp_xagg.readthedocs.io/)
[![Tests](https://github.com/rmcd-mscb/grd2shp_xagg/workflows/Tests/badge.svg)](https://github.com/rmcd-mscb/grd2shp_xagg/actions?workflow=Tests)
[![Codecov](https://codecov.io/gh/rmcd-mscb/grd2shp_xagg/branch/main/graph/badge.svg)](https://codecov.io/gh/rmcd-mscb/grd2shp_xagg)

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Features2

- TODO

## Requirements

- TODO

## Installation

You can install _Grd2Shp_Xagg_ via [pip](https://pip.pypa.io/) from
[PyPI](https://pypi.org/):

```console
$ pip install grd2shp_xagg
```

## Usage

Please see the [Command-line Reference](#usage) for details.

## Contributing

Contributions are very welcome. To learn more, see the [Contributor
Guide](CONTRIBUTING.rst).

## License

Distributed under the terms of the [MIT
license](https://opensource.org/licenses/MIT), _Grd2Shp_Xagg_ is free
and open source software.

## Issues

If you encounter any problems, please [file an
issue](https://github.com/rmcd-mscb/grd2shp_xagg/issues) along with a
detailed description.

### Develop

```console
$ conda env create -f environment.yml
$ conda develop -n {{cookiecutter.project_name}} src
$ conda activate {{cookiecutter.project_name}}
$ pip install -r requirements.dev
```

It is important to get [pre-commit]() enabled on the project, to ensure
that certain standards are always met on a git commit. With several of
these, it might fail if files are changed, but it will change them, and
trying the commit a second time will actually work.

### Git hook configuration

```console
$ pre-commit install --install-hooks
```

### Testing

[Nox]() is used for testing everything, with several sessions built-in.
To run the full suite of tests, simply use:

```console
$ nox
```

The different sessions are:

    - `pre-commit` -- validates that the [pre-commit]() checks all come back clean.
    - `safety` -- validates the [Safety]() of all production dependencies.
    - `mypy` -- validates the type-hints for the application using [mypy]().
    - `tests` -- runs all [pytest]() tests.
    - `typeguard` -- runs all [pytest]() tests, validates with [Typeguard]().
    - `xdoctest` -- runs any and all documentation examples with [xdoctest]().
    - `docs-build` -- builds the full set of generated API docs with [Sphinx]().

These can be run individually with the following command:

```console
$ nox -s <session>
```

Replace `<session>` with the name of the session give above, i.e.:

```console
$ nox -s mypy
```

You can also simply run [pytest]() tests, by using the command:

```console
$ pytest tests
```

### Dependencies

Production dependencies are duplicated, in both `requirements.txt` and
`environment.yml` due to how [conda]() does not work with the
`requirements.txt` file. It is necessary for both files to be updated as
dependencies are added.

Development dependencies are contained in `requirements.dev`.

### Version Management

    - Warning:
      Do not set the initial version to `0.0.0`. This will break, due to
      an unknown bug in [Bump2version]() at this time.

The projects made by this cookiecutter use [Bump2version]() for version
management. The default version that the project starts with is a
developmental version, `0.0.1-dev0`. In github, this should be
auto-incremented on each commit to the next dev build number. To manage
the version changes yourself, you can use the [Bump2version]() command:

```console
$ bump2version <part>
```

Where `<part>` is one of:

    - `major`
    - `minor`
    - `patch`
    - `build`

  <!-- end list -->

    - Note:
      This makes a `dev` version, which does not write a tag into git. It
      is just useful for development purposes and not the version that is
      recommended for a release version. The version string will be
      formatted as: `<major>.<minor>.<patch>-dev<build>`

To do a production release, use the command:

```console
$ bump2version --tag release
```

This will add a tag in the git repository noting the version.

    - Note:
      The version string for this will be: `<major>.<minor>.<patch>`

The template supports Python 3.8 and 3.9.

## Credits

This project was generated from
[@hillc-usgs](https://github.com/hillc-usgs)'s [Pygeoapi Plugin
Cookiecutter](https://code.usgs.gov/wma/nhgf/pygeoapi-plugin-cookiecutter)
template.
