Contributions
=============

Thank you for considering contributing to Plonk. Please read carefully the following guidelines on contributing to Plonk.

Code of conduct
---------------

We expect users and contributors to follow the [Python Community Code of Conduct](https://www.python.org/psf/codeofconduct/) in all modes of communication about Plonk. This means you are open, considerate, and respectful.

Getting started
---------------

If you want to contribute to Plonk you should fork the repository. You can then clone it to your local machine.

```bash
git clone https://github.com/your_user_name/plonk
```

Replace `your_user_name` with your GitHub user name.

There is a compiled Fortran component to Plonk which is derived from Splash. You must compile this before development. This requires a Fortran compiler, e.g. gfortran. We compile Splash into a shared object library, and then uses Cython to build a Python interface to that library. Then we use Conda to install Plonk in development mode. So, the steps are

1. Compile Splash into a library.
2. Build the Cython extension.
3. Install Plonk in Conda development mode.

There is a Makefile in the root directory of the repository to facilitate this.

```bash
make development
```

Then you can make changes to your local copy of Plonk and these changes will be reflected when you import Plonk and use it.

The Python interpreter must know where the Splash shared object library, `libsplash.so`, is at runtime. By default the Makefile installs it to `~/anaconda/lib`. If your Conda installation is different you can set the Makefile variable `INSTALL_DIR` as required

```bash
INSTALL_DIR=/path/to/conda/lib make development
```

You need to make sure the required dependencies are installed (via Conda). To satisfy these requirements there is a `environment.yml` file. You can set up a Conda environment for development and install Plonk in it. In the root directory of the repository, do the following

```bash
conda env create --file environment.yml
conda activate plonk-dev
```

and then follow the instructions above. To leave the development environment: `conda deactivate`.

If you make changes to Plonk that you would like to contribute you need to test that your code passes the test suite, and satisfies the code style requirements. To run the tests and check the code formatting do the following

```bash
make test
```

If this fails then there is either a test failure, or the code needs to be formatted in line with the chosen code style for Plonk. (See below.)

After you have committed and pushed your changes to your forked repository you
can issue a [pull request](https://github.com/dmentipl/plonk/pull/new/master).

Code style
----------

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for code style, and use [Black](https://github.com/python/black) and [isort](https://github.com/timothycrosley/isort) for auto-formatting.

> isort your python imports for you so you don't have to.

 Black is sponsered by the Python Software Foundation.

> Black is the uncompromising Python code formatter. By using it, you agree to cede control over minutiae of hand-formatting. In return, Black gives you speed, determinism, and freedom from pycodestyle nagging about formatting. You will save time and mental energy for more important matters.

To format your changes run the following from the main repository directory:

```bash
isort -rc
black --skip-string-normalization .
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

> ```
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

The documentation is not comprehensive.

- Add documentation for general usage.
- Add documentation for `plonk.Visualization`.
- Add documentation for `plonk.analysis`.

Documentation of use cases is also encouraged.

### Testing

We encourage contributions to the testing framework. To see where test coverage is lacking run `make test` and then open htmlcov/index.html in a web browser.

Coming up with ideas of what to test is useful.

### Features

Suggestions for new features include:

- additional analysis functions;
- a framework for modifying dump files;
- handling extra physics, such as magnetic fields, dust, and binary discs, etc.
