import numpy as np
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

from mzqc import MZQCFile as mzqc
import pyopenms as oms

"""
Calculate noise related QC metrics
"""


def getSN_medianmethod(spec: oms.MSSpectrum, norm: bool = True) -> float:
  """
    getSN_medianmethod get S/N from a spectrum via the median method

    As calculated in signal processing, signal and noise are descerned by the median
    intensity. The ratio is formed by the intensity sum in each category, scaled by
    the number of peaks recieved in each category.

    Parameters
    ----------
    spec : oms.MSSpectrum
        Spectrum to compute S/N from
    norm : bool, optional
        Scale by the number of peaks recieved in each category, by default True

    Returns
    -------
    float
        The ratio of signal-to-noise intensities
    """
  if spec.size() == 0:
    return 0.0

  median: float = 0.0
  maxi: float = 0.0
  spec.sortByIntensity(False)

  mar = np.array([s.getIntensity() for s in spec])
  median = np.median(mar)

  if not norm:
    return np.max(mar) / median

  sig = np.sum(mar[mar <= median]) / mar[mar <= median].size
  noi = np.sum(mar[mar > median]) / mar[mar > median].size
  # def sz():
  #     test = np.random.rand(30)
  #     median = np.median(test)
  #     sig = np.sum(test[test<=median])/test[test<=median].size

  # def ln():
  #     test = np.random.rand(30)
  #     median = np.median(test)
  #     sig = np.sum(test[test<=median])/len(test[test<=median])

  # from timeit import timeit
  # import numpy as np
  # timeit(sz, number=100000)
  # timeit(ln, number=100000)

  return sig / noi


def getSNMetrics(spectrum_acquisition_metrics_MS: mzqc.QualityMetric, ms_level: int) -> List[mzqc.QualityMetric]:
  """
    getSNMetrics collect S/N related QC metrics from a super metric collected in a first pass of the input mzML

    S/N from each spectrum are computed into 'spectrum acquisition metrics' for each MS level, from there S/N
    distribution values are computed.

    Parameters
    ----------
    spectrum_acquisition_metrics_MS : mzqc.QualityMetric
        QualityMetric object with the spectrum acquisition metrics
    ms_level : int
        The MS level to which the given spectrum acquisition metrics belong to

    Returns
    -------
    List[mzqc.QualityMetric]
        A list of new QualityMetric objects for mzQC deposition
    """
  metrics: List[mzqc.QualityMetric] = list()
  np_sn = np.array(spectrum_acquisition_metrics_MS.value['SN'])

  qs = np.quantile(np_sn, [.25, .5, .75])
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Signal-to-noise ratio Q1, Q2, Q3 for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=list(qs))
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Signal-to-noise ratio sigma for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=np.std(np_sn))
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Signal-to-noise ratio mean for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=np.mean(np_sn))
                 )

  low_out = qs[0] - (1.5 * (qs[2] - qs[0]))
  high_out = qs[2] + (1.5 * (qs[2] - qs[0]))
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Signal-to-noise ratio +/-1.5*IQR outlier for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=np.extract((np_sn < low_out) | (np_sn > high_out), np_sn))
                 )

  return metrics
