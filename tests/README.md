Running the tests
=================

From vscode with settings:
{
    "python.testing.unittestEnabled": true,
    "python.testing.unittestArgs": [
        "-v",
        "-s",
        "./tests",
        "-p",
        "test_*.py"
    ],
    "python.testing.pyTestEnabled": false,
    "python.testing.pyTestArgs": [
        "tests"
    ],
    "python.testing.nosetestsEnabled": false,
    "python.testing.autoTestDiscoverOnSaveEnabled": true,
    "python.testing.promptToConfigure": true,
}

From the command line (from the root folder of this repository):
python -m unittest discover -v -s ./tests -p test_*.py

Naming conventions
==================

test files: test_*.py
unittest class names: Test*
unittest class test methods: test_*
