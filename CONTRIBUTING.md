Contributions
=============

Thank you for considering contributing to Plonk. Please read carefully the following guidelines on contributing to Plonk.

Code of conduct
---------------

We expect users and contributors to follow the [Python Community Code of Conduct](https://www.python.org/psf/codeofconduct/) in all modes of communication about Plonk. This means you are open, considerate, and respectful.

Getting started
---------------

If you want to contribute to Plonk you should fork the repository. You can then clone it to your local machine, and use Conda to link to your local copy of the code.

```bash
git clone https://github.com/<user>/plonk
cd plonk && conda develop .
```

There is a compiled Fortran component to Plonk which is derived from Splash. You must compile this before development. This requires a Fortran compiler, e.g. gfortran. The following compiles Splash into a shared object library, and then uses Cython to build a Python interface to that library.

```bash
make install
python setup.py build_ext --inplace
```

This installs the libsplash shared object library at `~/anaconda/lib`. You will need to set `INSTALL_DIR` if your Python distribution is located elsewhere. (I.e. `make install INSTALL_DIR=your_python_lib_dir`.)

You need to make sure the required dependencies are installed (via Conda). To satisfy these requirements there is a `environment.yml` file. You can set up a Conda environment for development and install Plonk in it:

```bash
git clone https://github.com/<user>/plonk && cd plonk
conda env create --file environment.yml
conda activate plonk-dev
```

and then follow the instructions above. (To leave the development environment: `conda deactivate`.)

After you have committed and pushed your changes to your forked repository you
can issue a pull request: https://github.com/dmentipl/plonk/pull/new/master.

Code style
----------

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for code style, and use [Black](https://github.com/python/black) and [isort](https://github.com/timothycrosley/isort) for auto-formatting.

> isort your python imports for you so you don't have to.

 Black is sponsered by the Python Software Foundation.

> Black is the uncompromising Python code formatter. By using it, you agree to cede control over minutiae of hand-formatting. In return, Black gives you speed, determinism, and freedom from pycodestyle nagging about formatting. You will save time and mental energy for more important matters.

To format your changes run the following from the main repository directory:

```bash
isort plonk/**/*.py
black --skip-string-normalization plonk
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

The documentation is not very comprehensive.

- Add documentation for general usage.
- Add documentation for `plonk.Visualization`.
- Add documentation for `plonk.analysis`.

Documentation of use cases is also encouraged.

### Testing

The testing framework is not at all comprehensive. Contributions to the testing framework are strongly encouraged. For example, we require tests for

- reading and writing of Phantom dumps;
- calculating extra quantities;
- analysis functions;
- visualization functions.

Even just coming up with ideas of what to test is useful.

### Features

Suggestions for new features include:

- additional analysis functions;
- a framework for modifying dump files;
- handling extra physics, such as magnetic fields, dust, and binary discs, etc.
