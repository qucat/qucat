import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="qucat",
    version="0.1.1",
    author="Mario Gely",
    author_email="mario.f.gely@gmail.com",
    description="QUantum Circuit Analysis Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mgely/qucat",
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)