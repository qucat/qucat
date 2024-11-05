..  _contributing:

**********************
Contributing to QuCAT
**********************

In order to contribute
`create a pull request <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork>`_
from a fork of our 
`Github repository <https://github.com/qucat/qucat/>`_.
You can always :ref:`contact` us, 
especially if you want to bring big contributions to the project.

The project repository is organized as follows

- The source code is stored in the ``src`` folder, and most functions a user interacts with are in ``src/core.py``
- The functions and classes a user has access to when calling ``import qucat`` are defined in ``src/__init__.py``
- Unittesting is carried out automatically each time code is contributed through the testing scripts located in ``tests``. This is setup in the ``.github/workflows/tests.yml`` file
- The source code for the documentation is located in ``docs_src``. By executing the ``build_docs.py`` script, it becomes the content for the documentation website. This content is stored in the folder ``docs``
- One can use the cloned QuCAT library to test changes (rather than any pip-installed version) by editing and running the ``test.py`` script
- The ``master`` branch should reflect the latest pip-installable version of the software. Branches based on the ``master`` branch should be used for work in progress.

If you implement any new features in QuCAT, 
new function should contain 
a docstring and some comments, and
tests should be implemented to automatically verify your contribution.
Additional online documentation may be needed.
Tutorials to help in all these tasks can be found below

.. toctree::
   :maxdepth: 3

    Format for docstrings and comments <comments>
    How to write and run unittests <unittests>
    How to edit and build the documentation and website <online_doc>

