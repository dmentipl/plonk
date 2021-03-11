Contributions
=============

Thank you for considering contributing to Plonk. Please read the following guidelines on contributing to Plonk.

Code of conduct
---------------

We expect users and contributors to follow the [Python Community Code of Conduct](https://www.python.org/psf/codeofconduct/) in all modes of communication about Plonk. This means you are open, considerate, and respectful.

Getting started
---------------

If you want to contribute to Plonk you should first fork the repository. You can then clone it to your local machine.

```bash
# clone via HTTPS
git clone https://github.com/your_user_name/plonk

# or clone via SSH
$ git clone git@github.com:your_user_name/plonk
```

Replace `your_user_name` with your GitHub user name.

Set up an environment for Plonk development with Conda.

```bash
conda env create --file environment.yml
```

Use this environment for Plonk development.

```bash
conda activate plonk
```

Then you must create an "editable" install with pip.

```bash
cd plonk
python -m pip install -e .
```

Then you can make changes to your local copy of Plonk and these changes will be reflected when you import Plonk and use it. You can leave the development environment when done.

```bash
conda deactivate
```

If you make changes to Plonk that you would like to contribute, you need to test that your code passes the test suite, and satisfies the code style requirements. To run the tests do the following

```bash
python -m pytest
python -m mypy .
```

To check the code formatting do the following

```bash
python -m isort --check .
python -m black --check .
```

If any of these commands fail then there is either a test failure, or you need to reformat the code in line with the chosen code style for Plonk. (See below.)

After you have committed and pushed your changes to your forked repository you can issue a [pull request](https://github.com/dmentipl/plonk/pull/new/master).

To check the [code coverage](https://en.wikipedia.org/wiki/Code_coverage) do the following

```bash
python -m coverage run
python -m coverage html
```

Then open the report `htmlcov/index.html` in a web browser.

Code style
----------

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for code style, and use [Black](https://github.com/python/black) and [isort](https://github.com/timothycrosley/isort) for auto-formatting.

> isort your python imports for you so you don't have to.

 Black is sponsored by the Python Software Foundation.

> Black is the uncompromising Python code formatter. By using it, you agree to cede control over minutiae of hand-formatting. In return, Black gives you speed, determinism, and freedom from pycodestyle nagging about formatting. You will save time and mental energy for more important matters.

To format your changes run the following from the main repository directory:

```bash
python -m isort .
python -m black .
```

### Commit messages

When writing commit messages it's important to follow conventions to write clear and informative messages. We follow the advice contained in [The seven rules of a great Git commit message](https://chris.beams.io/posts/git-commit/#seven-rules):

> 1. Separate subject from body with a blank line
> 2. Limit the subject line to 50 characters
> 3. Capitalize the subject line
> 4. Do not end the subject line with a period
> 5. Use the imperative mood in the subject line
> 6. Wrap the body at 72 characters
> 7. Use the body to explain what and why vs. how

Here is an example from the [Git documentation](https://git-scm.com/book/ch5-2.html)

> ```git
> Short (50 chars or less) summary of changes
>
> More detailed explanatory text, if necessary.  Wrap it to about 72
> characters or so.  In some contexts, the first line is treated as the
> subject of an email and the rest of the text as the body.  The blank
> line separating the summary from the body is critical (unless you omit
> the body entirely); tools like rebase can get confused if you run the
> two together.
>
> Further paragraphs come after blank lines.
>
>   - Bullet points are okay, too
>
>   - Typically a hyphen or asterisk is used for the bullet, preceded by a
>     single space, with blank lines in between, but conventions vary here
> ```

Welcome contributions
---------------------

As well as adding new features, contributions to documentation and testing are welcome.

### Documentation

The documentation is not comprehensive. Documentation of use cases is encouraged.

### Testing

We welcome contributions to the testing framework. To see where test coverage is lacking run `python -m coverage run && python -m coverage html` and then open `htmlcov/index.html` in a web browser.

Coming up with ideas of what to test is useful.

### Features

Suggestions for new features include:

- [x] better support for physical units (with Pint);
- [ ] more analysis functions, e.g. for binary discs;
- [ ] a framework for modifying snapshot files;
- [ ] handling extra Phantom header quantities;
- [ ] out-of-core processing, e.g. using Dask or Vaex;
- [ ] extra visualization features, e.g. widgets in a Jupyter notebook with Bokeh;
- [ ] tracking particles through multiple snapshots;
- [ ] handling extra physics, such as magnetic fields.

New releases and PyPI and Conda packages
----------------------------------------

**Note: these instructions are for the Plonk maintainers (e.g. [@dmentipl](https://github.com/dmentipl)).**

### Make a release

First, increase the version number in `setup.cfg`, and update the `CHANGELOG.md` adding a heading like `[v0.3.1] - yyyy-mm-dd` under "Unreleased" heading. Then commit the change with a message like "Bump version to v0.3.1".

Then, make a new release on GitHub at <https://github.com/dmentipl/plonk/releases>. The title and tag should both be like "v0.3.1" which corresponds to the Plonk version number. Copy in the changes from the `CHANGELOG.md`. This creates a git tag for the commit, and generates a GitHub release with downloadable source as a tar.gz file.

### PyPI and pip

To generate a package installable from PyPI, first install twine (via the PyPI or Conda package). Then build a source and wheel distribution.

```bash
python setup.py sdist bdist_wheel
```

You can check that the build succeeded with the following.

```bash
python -m twine check dist/*
```

Then upload the package to the test PyPI to check that everything looks correct. You will be required to use your credentials.

```bash
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

If everything is good, upload to PyPI. The credentials are separate from the test PyPI.

```bash
python -m twine upload dist/*
```

### Conda

*Note: We use conda-forge to build the conda package.*

Clone [my fork](https://github.com/dmentipl/plonk-feedstock) of the Plonk feedstock. Modify the `meta.yml` file in two ways:

1. Update the version number.
2. Update the sha256 hash to correspond to the source version on [PyPI](https://pypi.org/project/plonk/).

Commit the change with a message like "Update to version 0.3.1". Then go to the GitHub page and generate a new [pull request](https://github.com/dmentipl/plonk-feedstock/pull/new/master). This will run several tests. If they pass, merge the pull request into the conda-forge/plonk-feedstock repository. Then a new conda package should soon be available on the [Anaconda cloud](https://anaconda.org/conda-forge/plonk).
