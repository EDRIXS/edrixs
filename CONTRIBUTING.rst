============
Contributing
============

Contributions are welcome and greatly appreciated. Bug reports, fixes,
features, tests, and documentation improvements all help the project.

We are currently updating parts of the EDRIXS infrastructure on the
``petsc_edrixs`` branch. Planning is tracked in the
`GitHub wiki <https://github.com/EDRIXS/edrixs/wiki>`_.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/EDRIXS/edrixs/issues.

If you are reporting a bug, please include:

* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

EDRIXS could always use more documentation, whether in the official
documentation, docstrings, examples, or tutorials.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at
https://github.com/EDRIXS/edrixs/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project.

Get Started!
------------

Ready to contribute? EDRIXS requires Python 3.10 or later. The compiled
extension also requires the native build dependencies described in the
`installation guide <docs/source/user/installation.rst>`_.

1. Fork the EDRIXS repository on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your-name/edrixs.git
    $ cd edrixs

3. Create a virtual environment and install EDRIXS with its development
   dependencies::

    $ python -m venv .venv
    $ source .venv/bin/activate
    $ python -m pip install --upgrade pip
    $ python -m pip install -r requirements-dev.txt
    $ python -m pip install -e .

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. Install the pre-commit hooks::

    $ pre-commit install

   The hooks run checks such as ``flake8`` when you commit. To check all files
   manually, run::

    $ pre-commit run --all-files

6. Run the test suite as described below, then commit and push your changes::

    $ git add .
    $ git commit -m "Your detailed description of your changes"
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Testing
-------

Run the complete test suite from the repository root::

    $ python -m pytest

During development, you can run a single test file or test by node ID::

    $ python -m pytest tests/test_solvers_dispatch.py
    $ python -m pytest tests/test_solvers_dispatch.py::test_build_op_rejects_unknown_backend

Tests that exercise integration between solver components use the
``integration`` marker. Run only those tests, or exclude them, with::

    $ python -m pytest -m integration
    $ python -m pytest -m "not integration"

Register any new custom pytest marker in ``pyproject.toml`` so that pytest can
validate it and does not emit ``PytestUnknownMarkWarning``.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests for behavior changes and bug fixes.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. Run the full test suite and pre-commit checks before submitting the pull
   request.
4. Make sure the continuous-integration checks pass on all supported Python
   versions.
