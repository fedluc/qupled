How to contribute
=================

The following guidelines explain how to contribute to both the codebase and the documentation. To properly set up the
tools mentioned below, start by installing the required Python packages with the following command:

.. code-block:: console

   pip install -r dev/requirements.txt

Developing
----------

Qupled is an open-source project, and contributions from external developers are always welcome.
If you'd like to add or test new features, you can build the code by running:

.. code-block:: console

   ./devtool build

This command, executed from the root directory of the project, will generate the build and place it in the ``dist``
directory as a Python wheel file. To verify that everything is working correctly, run:

.. code-block:: console

   ./devtool test

Tests are also executed automatically whenever new code is pushed to the remote repository. If all the tests
pass correctly you can then proceed to install your new version of qupled with:

.. code-block:: console

   pip install --force-reinstall dist/qupled*.whl

Formatting
----------

To maintain consistent code formatting across the C++ and Python codebases, we use
`clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ and
`black <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html>`_.
The formatting is automatically checked every time new code is pushed to the repository.
To manually ensure the correct formatting is applied, run:

.. code-block:: console

   ./devtool format

Documentation
-------------

The documentation is stored in the ``docs`` directory, and changes can be made by editing the ``.rst`` files within it.
Once you've made your changes, you can verify and build the documentation using:

.. code-block:: console

   ./devtool docs

The generated output can be viewed by opening ``docs/_build/index.html`` in your browser.
