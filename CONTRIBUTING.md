Contributions
=============

Thank you for considering contributing to Plonk. Please read the following guidelines on contributing to Plonk.

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

Set up an environment for Plonk development with Conda.

```bash
conda env create --file environment.yml
conda develop --name plonk-dev .
```

Use this environment for Plonk development.

```bash
conda activate plonk-dev
```

Then you can make changes to your local copy of Plonk and these changes will be reflected when you import Plonk and use it. You can leave the development environment when done.

```bash
conda deactivate
```

If you make changes to Plonk that you would like to contribute, you need to test that your code passes the test suite, and satisfies the code style requirements. To run the tests and check the code formatting do the following

```bash
python -m coverage run -m pytest && coverage html
isort --skip plonk/__init__.py --check-only -rc
black --check --skip-string-normalization plonk tests
mypy --no-strict-optional --ignore-missing-imports plonk tests
```

If any of these commands fail then there is either a test failure, or you need to reformat the code in line with the chosen code style for Plonk. (See below.)

After you have committed and pushed your changes to your forked repository you
can issue a [pull request](https://github.com/dmentipl/plonk/pull/new/master).

Code style
----------

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for code style, and use [Black](https://github.com/python/black) and [isort](https://github.com/timothycrosley/isort) for auto-formatting.

> isort your python imports for you so you don't have to.

 Black is sponsored by the Python Software Foundation.

> Black is the uncompromising Python code formatter. By using it, you agree to cede control over minutiae of hand-formatting. In return, Black gives you speed, determinism, and freedom from pycodestyle nagging about formatting. You will save time and mental energy for more important matters.

To format your changes run the following from the main repository directory:

```bash
isort --skip plonk/__init__.py -rc
black --skip-string-normalization plonk tests
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

We welcome contributions to the testing framework. To see where test coverage is lacking run `python -m coverage run -m pytest && coverage html` and then open `htmlcov/index.html` in a web browser.

Coming up with ideas of what to test is useful.

### Features

Suggestions for new features include:

- better support for physical units (with Pint);
- more analysis functions, e.g. for binary discs;
- a framework for modifying snapshot files;
- handling extra Phantom header quantities;
- out-of-core processing, e.g. using Dask or Vaex;
- extra visualization features, e.g. widgets in a Jupyter notebook with Bokeh;
- tracking particles through multiple snapshots;
- handling extra physics, such as magnetic fields.
