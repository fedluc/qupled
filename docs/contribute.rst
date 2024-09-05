How to contribute
=================

If you want to make contributions to the code make sure to read the following guidelines concerning the
formatting and documentation.

Formatting
----------

We use `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ and `black <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html>`_ to enforce a uniform formatting style for the C++ and Python code-base. The code formatting is checked automatically every time a pull request or a push is made. To ensure that the correct formatting is applied run ``cmake --build . --target format`` from the build directory.

Documentation
-------------

All the documentation is stored in the ``docs`` folder. Changes can be made by editing the ``.rst`` files in ``docs``. After the changes are applied the documentation should be built in order to check that it has the correct format. The `sphinx <https://www.sphinx-doc.org/en/master/>`_ and `sphinx_rdt_theme <https://pypi.org/project/sphinx-rtd-theme/>`_ python packages must be installed in order to build the documentation correctly. The documentation is built by navigating to the ``docs`` folder and running

.. code-block:: console

   sphinx-build -b html . _build

The result of the build can be inspected by opening ``docs/_build/index.html``.
