How to contribute
=================

The following guidelines explain how to contribute to both the codebase and the documentation. 

Setup the development environment
---------------------------------

There are two options for setting up your development environment:

1. **The Easy Option: Use the Development Container**

   This project includes a pre-configured development container to simplify the setup process. 
   The development container provides all necessary tools, including an up-to-date version of Git, 
   Python, and other dependencies.

   To use the development container, ensure you have Docker and Visual Studio Code installed. 
   Then, open the project in Visual Studio Code and follow these steps:

   - Install the `Remote - Containers` extension if you haven't already.
   - When prompted, reopen the project in the container.

   Once the container is running, all tools and dependencies will be available in the environment. 
   You can start developing immediately without additional setup.

2. **The DIY Option: Manual Setup**

   If you prefer to set up the environment manually, install `uv <https://docs.astral.sh/uv/>`_
   first and bootstrap ``foga`` with:

   .. code-block:: console

      uv sync --group dev --no-install-project

   Then use ``foga`` to install the required system dependencies:

   .. code-block:: console

      .venv/bin/foga install --target apt-get   # Ubuntu or Debian-based systems
      .venv/bin/foga install --target brew      # macOS

   and sync the full Python development environment:

   .. code-block:: console

      .venv/bin/foga install --target dev-env

   This keeps ``uv`` as the bootstrap layer and lets ``foga`` own the project-specific
   environment setup. After that, use ``foga`` for the repository workflows:

   .. code-block:: console

      foga build
      foga test --runner unit

   Additionally, ensure that you have all the necessary
   :ref:`external dependencies <external_dependencies>` installed.

Formatting
----------

To maintain consistent code formatting across the C++ and Python codebases, we use
`clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ and
`black <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html>`_.
The formatting is automatically checked every time new code is pushed to the repository.
To manually ensure the correct formatting is applied, run:

.. code-block:: console

   foga format

Testing
-------

Before opening a pull request, run the test suite for the area you changed.

For Python tests:

.. code-block:: console

   foga test --runner unit
   foga test --runner native
   foga test --runner integration

For C++ tests:

.. code-block:: console

   foga test cpp --runner cpp

MPI is disabled by default for ``cpp`` tests; use the ``mpi`` profile to enable it.

.. code-block:: console

   foga test cpp --runner cpp --profile mpi

To run everything in one shot:

.. code-block:: console

   foga test

Documentation
-------------

The documentation is stored in the ``docs`` directory, and changes can be made by editing the ``.rst`` files within it.
Once you've made your changes, you can verify and build the documentation using:

.. code-block:: console

   foga docs

The generated output can be viewed by opening ``docs/_build/index.html`` in your browser.
