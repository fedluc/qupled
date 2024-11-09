How to contribute
=================

If you want to make contributions to the code make sure to read the following guidelines concerning the
formatting and documentation.

Formatting
----------

We use `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ and `black <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html>`_ to enforce a uniform formatting style for the C++ and Python code-base. The code formatting is checked automatically every time a pull request or a push is made. To ensure that the correct formatting is applied run ``cmake --build . --target format`` from the build directory.

Documentation
-------------

The documentation is located in the ``docs`` folder, with changes made by editing the ``.rst`` files within it. Once modifications are complete, the documentation should be rebuilt to ensure it is correctly formatted. To build the documentation, you’ll need the Python packages `sphinx <https://www.sphinx-doc.org/en/master/>`_ and `sphinx_rdt_theme <https://pypi.org/project/sphinx-rtd-theme/>`_. From the build directory, use the command ``cmake --build . --target docs`` to build the documentation, which will also compile the qupled library if it hasn't been built yet. The build output can be viewed by opening ``docs/index.html`` in the build directory.
