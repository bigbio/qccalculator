from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='qccalculator',
    version='0.0.2',
    packages=find_packages(),
    url='https://gitlab.ebi.ac.uk/walzer/qccalculator',
    description='Development for python driven QC calculation (qccalculator) and presentation (QSpector)',
    long_description=long_description,
    install_requires=[
        "biopython",
        "click",
        "mzqc-pylib",
        "pandas",
        "plotly-express",
        "pronto",
        "pyopenms",
        "requests",
        "toposort"
        ],
    include_package_data=True,
    #scripts=['bin/qspector/qspector.py'] ,
    entry_points={
        'console_scripts': ['qccalculator = qccalculator.cli:start']
      }
)
