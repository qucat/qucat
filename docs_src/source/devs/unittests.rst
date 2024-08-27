.. _unittesting:

=========================
Unit testing
=========================

Tests are setup using pythons default testing framework 
`unittest <https://docs.python.org/2/library/unittest.html>`_.
These tests are run each time new code is pushed to our Github repository
using `GitHub Actions <https://docs.github.com/en/actions/>`_.
They can also be edited
and run on your local machine as explained below.

In these explanations, we assume you have 
cloned a version of the
`qucat repository <https://github.com/qucat/qucat/>`_
to your computer.

Running tests
=============

Tests are located in the the ``tests`` directory. These can be run from the command line 
after navigating to the root directory of the qucat repository with the command

``python -m unittest discover -v -s ./tests -p test_*.py``

Each test file is dedicated to a specific part of the library, 
and is generally paired with a specific file of the source code, 
for example the file ``test/test_core.py`` tests the functionalities of the 
``src/test_core.py`` file.

To run only run one of the test files, ``test_core.py`` for example, one can also run 

``python tests/test_core.py``

Executing only a single test from a file is also possible
by running for example

``python tests/test_core.py SeriesRLC.test_frequency``

.. note:: Some GUI tests have failed for a mysterious reason on certain
    operating systems, even thought the functionality they test works 
    when tested manually.

Writing a basic test
====================

Let us assume that a developer wishes to contribute a new module ``src/new.py``
which contains a function ``f`` which adds two numbers

.. literalinclude:: example_minimal_new_module.txt
   :language: python

One should then create a ``test/test_new.py`` file which test if the 
function works as expected.
This is one way of doing it

.. literalinclude:: example_minimal_test.txt
   :language: python

For more details on the unittest framework, check
`this link <https://docs.python.org/2/library/unittest.html>`_.

Writing a test for a GUI functionality
======================================

Testing the GUI functionalities *i.e.* **what 
happens if I click here?** is also possible in an automatic 
way by simulating mouse motion, clicks or keystrokes events.
These test are located in ``test/test_gui.py``

Creating these tests is also done through the GUI, and we present
here a short tutorial on how to do so.

Let's assume we want to test moving a resistor by clicking, dragging
then dropping.

The correct location to add this test is in the 
``test/test_gui.py`` file, with the following code

.. literalinclude:: example_gui_test.txt
   :language: python

The **first time** we run this file,
an empty GUI will open.
**Step 1: ** here we set the initial configuration for the
test, in this example we create and place a resistor, then close the GUI.
This circuit will be saved in the file 
``tests/gui_testing_files/test_moving_resistor/initial_netlist.txt``
**Step 2: **  a second GUI will open in the initial configuration, and we now 
perform the task we are testing. 
Our actions will be recorded in the file
``tests/gui_testing_files/test_moving_resistor/events.txt`` 
later repeated in an automatic way when running this test. In this 
example, we drag and drop the resistor to another location, then
close the GUI.
Upon closing the GUI, the final configuration of the circuit
is recorded in the file
``tests/gui_testing_files/test_moving_resistor/final_netlist``

The **subsequent times** we run this test no human intervention
is necessary.
A GUI will be opened in the initial configuration, 
a sequence of events will be triggered, 
and the final configuration of the circuit will be stored in the file
``tests/gui_testing_files/test_moving_resistor/final_after_events_netlist.txt``
and compared to the configuration created manually
``tests/gui_testing_files/test_moving_resistor/final_netlist``.
If the files are identical, the test passes.

The function ``AutomaticTesting.launch_gui_testing`` which implements
these functionalities takes two arguments.
By setting ``force_build = True``, the test is run as if for the first 
time. 
By setting ``run_slower = True``, the test is run very slowly
so that one can observe what the automatically generated sequence 
of events is doing to the GUI. 
This enables easy debugging of the test.

Writing a test based on a tutorial
======================================

All tutorials featured on the website are 
jupyter notebooks
located in the folder ``docs_src/source/tutorials``
and are
also run as tests.
These tests are located in ``test/test_tutorials.py``.

If one adds a new notebook, for example ``new_tutorial.ipynb``
in the tutorials folder, one can test it by adding the following code 
to the ``test/test_tutorials.py`` script

.. literalinclude:: example_tutorial_test.txt
   :language: python