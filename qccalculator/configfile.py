# NOTE: do not change any values in here!

config_content = """[dalton.tolerances]
FTMS = 0.05
ITMS = 0.5
QTOF = 0.1
default = 0.5

[ppm.tolerances]
FTMS = 0.05
ITMS = 0.5
QTOF = 0.1
default = 0.5

[fdr.cutoff]
default = 0.05
strict = 0.01

[onto.urls]
psi-ms = https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
psi-qc = https://raw.githubusercontent.com/HUPO-PSI/mzQC/master/cv/qc-cv.obo
units = http://purl.obolibrary.org/obo/uo.owl
pride = http://purl.obolibrary.org/obo/pride_cv.obo
"""