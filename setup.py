from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='qccalculator',
    version='0.0.3',
    packages=find_packages(),
    url='https://gitlab.ebi.ac.uk/walzer/qccalculator',
    description='Development for python driven QC calculation (qccalculator)',
    long_description=long_description,
    install_requires=[
        "mzqc-pylib>=1.0.0 @ git+https://github.com/bigbio/mzqc-pylib.git@v1.0.0",
        "biopython",
        "click",
        "pandas",
        "plotly-express",
        "pronto",
        "pyopenms",
        "requests",
        "toposort",
        "plotly"
        ],
    include_package_data=True,
    #scripts=['bin/qspector/qspector.py'] ,
    entry_points={
        'console_scripts': ['qccalculator = qccalculator.cli:cli', 'mzqcfileinfo = qccalculator.fileinfo:mzqcfileinfo']
      }
)
