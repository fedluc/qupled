In order to build the sphinx documentation first install the following packages:

- sphinx (``pip install sphinx``)
- sphinx_rtd_theme (``pip install sphinx_rtd_theme``)

Once changes have been made in the .rst files in order to see them in the
documentation go to qupled/docs directory and run:

sphinx-build -b html . _build

In the qupled/docs/_build directory open the index.html file in a browser.