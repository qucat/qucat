import setuptools

setuptools.setup(
    long_description='QUCAT: QUantum Circuit Analyzer Tool.\n\nSee https://qucat.org/ for installation, documentation, tutorials and more.',
    long_description_content_type="text/markdown",
    package_dir={'qucat': 'src'},
    packages= ['qucat'],
    include_package_data=True,
    package_data={'qucat': ['.graphics/*']},
)
