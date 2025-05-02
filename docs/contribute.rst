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

   If you prefer to set up the environment manually, start by installing the required Python 
   packages with the following command:

   .. code-block:: console

      pip install -r dev/requirements.txt

   Additionally, ensure that you have all the necessary :ref:`external dependenceis <external_dependencies>` 
   installed.

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
