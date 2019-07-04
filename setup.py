import setuptools
from distutils.dir_util import remove_tree

setuptools.setup(
    name="qucat",
    version="0.7.2",
    author="Mario Gely",
    author_email="mario.f.gely@gmail.com",
    description="QUantum Circuit Analysis Tool",
    long_description='QUCAT: QUantum Circuit Analyzer Tool.\n\nSee https://qucat.org/ for installation, documentation, tutorials and more.',
    long_description_content_type="text/markdown",
    url="https://qucat.org",
    package_dir={'qucat': 'src'},
    packages= ['qucat'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>3.0.0',
    include_package_data=True,
    package_data={'qucat': ['.graphics/*']},
)