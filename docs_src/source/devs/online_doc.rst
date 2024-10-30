========================================
Writing documentation or website content
========================================

We assume you have 
cloned a version of the
`qucat repository <https://github.com/qucat/qucat/>`_
to your machine.

The source code for the documentation is located in the 
``/docs_src/source`` folder.
The built documentation is stored locally in the ``/docs/`` folder. 
Through GitHub Pages, the website (qucat.org) always
matches the contents of the ``/(root)/`` folder of the 
``gh-pages`` branch.
The website can be previewed locally by opening the ``/docs/index.html`` 
file with a web browser.

In order to edit the documentation/website one should thus follow 
these steps

- edit the reStructuredText files located in ``/docs_src/source``
- build the documentation by running the ``build_docs.py`` script
- check the changes locally by opening the ``/docs/index.html`` file with a web browser
- push the code to the GitHub repository. Note that the ``/docs/`` folder is not pushed to the GitHub repository
- once the code is pulled into the master branch, a GitHub Action will build the website and push the new ``/docs/`` folder to the root of the ``gh-pages`` branch, which will automatically update the website.
Do not make changes to the ``gh-pages`` branch directly, as it is automatically updated by the GitHub Action.

Building the documentation,
and editing the tutorials or the documentation of the functions
requires specific instruction given below.

Building the documentation
--------------------------

The documentation can be built by running the ``build_docs.py`` script.

Some details about the build process: 
we use Sphinx to build the html content for the website
from the reStructuredText (.rst) files located in ``/docs_src/source``.
This build process is configured through the ``/docs_src/source/conf.py``.
Information on sphinx, on reStructuredText and the configuration file 
can be found online, 
notably `here <http://www.sphinx-doc.org/en/master/>`_.

.. note:: Building the docs requires `pandoc <https://pandoc.org/installing.html>`_, and the pip-installable packages sphinx, nbsphinx, recommonmark, sphinx_rtd_theme

Editing class or function documentation
--------------------------------------------

One advantage of using Sphinx, is that it can 
extract docstrings from the source-code
(located in the ``src`` folder) and turn
it into elegant html documentation of 
a class or function.

As an example for classes, inserting the following code
in the reStructuredText file of one of the documentation pages

.. literalinclude:: example_class_autostring.txt
   :language: python

will produce the content of the page 
:ref:`capacitor`.

Here ``:members:`` and ``:inherited-members:`` indicates that
we also want to display information about the class 
methods and inherited methods.
One can similarly include documentation for 
functions using ``.. automodule::``, more 
details on the topic can be found 
`here <http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_.

Since the automatic documentation 
has already been set up for all QuCAT classes,
**any public modules added to existing classes
by a contributor, and correctly commented, will thus be automatically
included in the documentation website after their pull request
is merged into the master branch**.
Details about how to write good docstrings 
which lead to 
such documentation can be found at 
:ref:`commenting`.

Note that one can also create a private module, 
which will not be included 
in the documentation,
by prepending the name
of the function by an underscore.

Editing the tutorials
--------------------------------------------

The tutorials featured on the website correspond
to jupyter notebooks located in the ``docs_src/source/tutorials``
folder.

In order to **edit** a tutorial

- Open a tutorial in ``docs_src/source/tutorials``, edit and save it

In order to **add a new** tutorial

- Add your jupyter notebook to the ``docs_src/source/tutorials`` folder.
- Edit the  ``docs_src/source/tutorials/index.rst`` file by adding the name of your notebook to the existing list

In **both cases**, subsequently

- Build the documentation by running the ``build_docs.py`` script
- Preview your changes by opening the ``/docs/index.html`` file with a web browser
- Make a pull request: your changes will be visible on the website when merged to the master branch

If you are adding a new tutorial, you may want to include it in the 
set of tests QuCAT performs each time a change is made to the code, see 
:ref:`unittesting`.
