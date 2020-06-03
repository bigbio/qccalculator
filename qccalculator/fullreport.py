from pyopenms import MSExperiment
from qccalculator.basicqc import getBasicQuality


class FullReport:

  def __init__(self, id_file: str ):
    # The spectra file is a dictionary with the file name -> OpenMS experiment
    self._spectra_files = {}
    self._id_file = id_file
    self._basicqc = {}

  def add_openms_experiment(self, spectra_file: str, msexperiment: MSExperiment):
    self._spectra_files[spectra_file] = msexperiment

  def basicqc(self):
    for file_key in self._spectra_files:
      basic = getBasicQuality(self._spectra_files[file_key], verbose=True)
      self._basicqc[file_key] = basic
