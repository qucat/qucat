#!/bin/bash

# Help function to publish to (test)PyPI
# Usage: ./publish_to_PyPI.sh [pypi|testpypi]
# If no argument is provided, the default is testPyPI
# The version number is read from the pyproject.toml file.
# The version number is prefixed with 'testv' for the testPyPI repository
# and 'v' for the (production) PyPI repository.

# Default repository
repository="testpypi"
# Set a prefix for the version
version_prefix="testv"

# Check if the repository is provided
if [ "$1" ]
then
    # check that the repository is either testpypi or pypi
    if [ "$1" = "pypi" ]
    then
        version_prefix="v"
    elif [ "$1" != "testpypi" ]
    then
        echo "Invalid repository. Please provide either 'testpypi' or 'pypi'"
        exit 1
    fi
fi

echo "Using repository: $repository"
git tag $version_prefix$(grep version pyproject.toml | awk -F'"' '$0=$2')
git push origin --tags
