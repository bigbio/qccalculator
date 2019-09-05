from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='QCCalculator',
    version='0.0.1',
    packages=find_packages(),
    url='https://gitlab.ebi.ac.uk/walzer/qccalculator',
    description='Development for python driven QC calculation (QCCalculator) and presentation (QSpector)',
    long_description=long_description,
    install_requires=[
        "rpy2",
        "MZQC",
        ],
    include_package_data=True,
    scripts=['QCCalculator.py'] ,
)