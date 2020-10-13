import datetime
import itertools
import sys
from collections import defaultdict
from os.path import basename
from typing import Any, Dict, List

import numpy as np
import pyopenms as oms
from pyopenms import ActivationMethod
from mzqc import MZQCFile as mzqc

from qccalculator import utils
from qccalculator import noiseqc

"""
Basic set of methods to get quality metric calculations from peak and identification files under way
"""


def getBasicQuality(exp: oms.MSExperiment, verbose: bool = False) -> mzqc.RunQuality:
  """
    getBasicQuality calculates the basic QualityMetrics from a mass spectrometry peak file and creates the related RunQuality object.

    Calculated basic QC metrics and proto-metrics necessary to calculate more elaborate QC metrics with additional data (e.g. ID).

    Parameters
    ----------
    exp : oms.MSExperiment
        The mass spectrometry peak file to calculate metrics from
    verbose : bool, optional
        switches on verbose logging, by default False

    Returns
    -------
    mzqc.RunQuality
        A RunQuality object containing the list of metrics calculated and metadata collected, ready for integration into a mzQC file object.
    """
  metrics: List[mzqc.QualityMetric] = list()
  if exp.getExperimentalSettings().getSourceFiles():
    parent_base_name: str = basename(exp.getExperimentalSettings().getSourceFiles()[0].getNameOfFile())
    parent_chksm: str = exp.getExperimentalSettings().getSourceFiles()[0].getChecksum()
    parent_chksm_type: str = exp.getExperimentalSettings().getSourceFiles()[0].getChecksumType()

  instr_srl: str = exp.getInstrument().getMetaValue('instrument serial number') \
    if exp.getInstrument().metaValueExists('instrument serial number') else 'unknown'  # MS:1000529 in mzML

  input_loc: str = exp.getExperimentalSettings().getLoadedFilePath()
  base_name: str = basename(input_loc)
  chksm: str = utils.sha256fromfile(exp.getExperimentalSettings().getLoadedFilePath())
  cmpltn: str = exp.getDateTime().get()
  # strt:datetime.datetime = datetime.datetime.strptime(cmpltn, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(seconds=exp.getChromatograms()[0][exp.getChromatograms()[0].size()-1].getRT()*60)

  meta: mzqc.MetaDataParameters = mzqc.MetaDataParameters(
    inputFiles=[
      mzqc.InputFile(name=base_name, location=input_loc,
                     fileFormat=mzqc.CvParameter("MS", "MS:1000584", "mzML format"),
                     fileProperties=[
                       mzqc.CvParameter(cvRef="MS",
                                        accession="MS:1000747",
                                        name="completion time",
                                        value=cmpltn
                                        ),
                       mzqc.CvParameter(cvRef="MS",
                                        accession="MS:1000569",
                                        name="SHA-256",
                                        value=chksm
                                        ),
                       mzqc.CvParameter(cvRef="MS",
                                        accession="MS:1000031",
                                        name="instrument model",
                                        value=exp.getInstrument().getName()
                                        ),
                       mzqc.CvParameter(cvRef="MS",
                                        accession="MS:1000529",
                                        name="instrument serial number",
                                        value=instr_srl
                                        )
                       # TODO integrate parent location and checksum
                       # id: MS:1002846 (Associated raw file URI) N.B. definition is PRIDE specific - WTF
                       # fitting checksum cv missing
                     ]
                     )
    ],
    analysisSoftware=[
      mzqc.AnalysisSoftware(cvRef="MS", accession="MS:1000752", name="TOPP software", version=oms.__version__,
                            uri="openms.de")
    ]
  )

  # this is mighty important to sort by RT
  exp.sortSpectra()

  min_mz: float = sys.maxsize
  max_mz: float = 0
  mslevelcounts: Dict[int, int] = defaultdict(int)

  spectrum_acquisition_metrics_MS1: Dict[str, List[Any]] = defaultdict(list)
  spectrum_acquisition_metrics_MS2: Dict[str, List[Any]] = defaultdict(list)
  spectrum_topn: Dict[str, List[Any]] = defaultdict(list)
  tandem_spectrum_metrics_MS2: Dict[str, List[Any]] = defaultdict(list)
  trap_metrics_MS1: Dict[str, List[Any]] = defaultdict(list)
  trap_metrics_MS2: Dict[str, List[Any]] = defaultdict(list)
  isolation_window_metrics: Dict[str, List[Any]] = defaultdict(list)
  tic_tab: Dict[str, List[Any]] = defaultdict(list)

  # ActivationMethod look-up dict
  ams = {getattr(ActivationMethod, i): i for i in dir(ActivationMethod) if type(getattr(ActivationMethod, i)) == int}

  intens_sum: np.float = 0
  last_surveyscan_index: int = 0
  last_surveyscan_intensity = -1
  last_surveyscan_max = -1

  for spin, spec in enumerate(exp):
    mslevelcounts[spec.getMSLevel()] += 1

    iontraptime = utils.getTrapTime(spec)
    intens_max = spec.get_peaks()[1].max()
    intens_min = spec.get_peaks()[1].min()
    intens_sum = spec.get_peaks()[1].sum()

    if spec.getMSLevel() == 1:
      last_surveyscan_index = spin
      last_surveyscan_intensity = intens_sum
      last_surveyscan_max = intens_max

      spectrum_acquisition_metrics_MS1['RT'].append(spec.getRT())
      spectrum_acquisition_metrics_MS1['SN'].append(noiseqc.getSN_medianmethod(spec))
      spectrum_acquisition_metrics_MS1['peakcount'].append(spec.size())
      spectrum_acquisition_metrics_MS1['int'].append(intens_sum.item())  # .item() for dtype to pytype

      trap_metrics_MS1['RT'].append(spec.getRT())
      trap_metrics_MS1['traptime'].append(iontraptime)

      tic_tab['RT'].append(spec.getRT())
      tic_tab['int'].append(intens_sum)

    if spec.getMSLevel() == 2:
      if spec.getPrecursors()[0].getMZ() < min_mz:
        min_mz = spec.getPrecursors()[0].getMZ()
      if spec.getPrecursors()[0].getMZ() > max_mz:
        max_mz = spec.getPrecursors()[0].getMZ()

      spectrum_acquisition_metrics_MS2['RT'].append(spec.getRT())
      spectrum_acquisition_metrics_MS2['SN'].append(noiseqc.getSN_medianmethod(spec))
      spectrum_acquisition_metrics_MS2['peakcount'].append(spec.size())
      spectrum_acquisition_metrics_MS2['int'].append(intens_sum.item())  # .item() for dtype to pytype
      spectrum_acquisition_metrics_MS2['native_id'].append(utils.spec_native_id(spec))

      rank = spin - last_surveyscan_index
      spectrum_acquisition_metrics_MS2['rank'].append(rank)

      trap_metrics_MS2['RT'].append(spec.getRT())
      trap_metrics_MS2['traptime'].append(iontraptime)
      trap_metrics_MS2['activation_method'].append(
        ams.get(next(iter(spec.getPrecursors()[0].getActivationMethods()), None), 'unknown'))
      trap_metrics_MS2['activation_energy'].append(spec.getPrecursors()[0].getMetaValue('collision energy') if \
                                                     spec.getPrecursors()[0].metaValueExists(
                                                       'collision energy') else -1)

      precursor_index = \
        np.searchsorted(exp[last_surveyscan_index].get_peaks()[0], [exp[spin].getPrecursors()[0].getMZ()])[0]
      if precursor_index != np.array(exp[last_surveyscan_index].get_peaks()).shape[1]:
        precursor_err = spec.getPrecursors()[0].getMZ() - \
                        np.array(exp[last_surveyscan_index].get_peaks())[:, precursor_index][0]
        precursor_int = np.array(exp[last_surveyscan_index].get_peaks())[:, precursor_index][1]
      else:
        precursor_err = np.nan
        precursor_int = np.nan

      tandem_spectrum_metrics_MS2['RT'].append(spec.getRT())
      tandem_spectrum_metrics_MS2['precursor_intensity'].append(
        precursor_int)  # TODO different from mzid->mzml getPrecursors[0].getIntensity() ? YES, latter one usually zero
      tandem_spectrum_metrics_MS2['precursor_error'].append(precursor_err)
      tandem_spectrum_metrics_MS2['precursor_mz'].append(spec.getPrecursors()[0].getMZ())
      tandem_spectrum_metrics_MS2['precursor_c'].append(spec.getPrecursors()[0].getCharge())

      if last_surveyscan_intensity != -1:
        tandem_spectrum_metrics_MS2['surveyscan_intensity_sum'].append(last_surveyscan_intensity)
      if last_surveyscan_max != -1:
        tandem_spectrum_metrics_MS2['surveyscan_intensity_max'].append(last_surveyscan_max)

      isolation_window_metrics['RT'].append(spec.getRT())
      isolation_window_metrics['isolation_target'].append(spec.getPrecursors()[
                                                            0].getMZ())  # https://github.com/OpenMS/OpenMS/blob/d17cc251fd0c4068eb253b03c9fb107897771fdc/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L1992
      isolation_window_metrics['isolation_lower'].append(spec.getPrecursors()[0].getIsolationWindowLowerOffset())
      isolation_window_metrics['isolation_upper'].append(spec.getPrecursors()[0].getIsolationWindowUpperOffset())
      lower = spec.getPrecursors()[0].getMZ() - spec.getPrecursors()[0].getIsolationWindowLowerOffset()
      upper = spec.getPrecursors()[0].getMZ() + spec.getPrecursors()[0].getIsolationWindowUpperOffset()

      s = np.array([(i.getMZ(), i.getIntensity()) for i in exp[last_surveyscan_index]], ndmin=2)
      s = s[np.where(np.logical_and(s[:, 0] >= lower, s[:, 0] <= upper))[0]]
      isolation_window_metrics['peaks_in_window'].append(np.shape(s)[0])

      int_sort_desc = np.flip(np.argsort(s[:, 1]))
      if np.shape(s)[0] > 1:
        isolation_window_metrics['int_ratio_ranked_peaks_in_window'].append(
          s[int_sort_desc][:-1, 1] / s[int_sort_desc][1:, 1][0])  # intensity ratio between top1&2, 2&3, ...
      else:
        isolation_window_metrics['int_ratio_ranked_peaks_in_window'].append(0)  # bigger is better, though best is 0

      isolation_window_metrics['summed_window_intensity'].append(np.sum(s[int_sort_desc][:, 1]))
      isolation_window_metrics['isolation_target_intensity'].append(spec.getPrecursors()[0].getIntensity())

      # TODO this needs to go outside
      tol = 0.5
      if spec.metaValueExists('filter string'):
        if 'FTMS' in spec.getMetaValue('filter string'):
          tol = 0.05
        elif 'ITMS' in spec.getMetaValue('filter string'):
          tol = 0.5
        elif 'QTOF' in spec.getMetaValue('filter string'):  # TOFMS, SQMS, TQMS, SectorMS
          tol = 0.1

      # ms2 peaks directly from isolation window?
      unfragmented = np.any([np.isclose(i[0], [x.getMZ() for x in spec], atol=tol) for i in s])
      isolation_window_metrics['peaks_in_window_in_ms2'].append(str(unfragmented))

  ## Spectra detail numbers
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Spectrum acquisition metric values - MS1",
                       value=spectrum_acquisition_metrics_MS1)
  )
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Spectrum acquisition metric values - MS2",
                       value=spectrum_acquisition_metrics_MS2)
  )
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Spectra topn ranks",
                       value=spectrum_topn)
  )
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Tandem spectrum metric values - MS2",
                       value=tandem_spectrum_metrics_MS2)
  )
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Trap metric values - MS1",
                       value=trap_metrics_MS1)
  )
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Trap metric values - MS2",
                       value=trap_metrics_MS2)
  )
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="isolation window metrics",
                       value=isolation_window_metrics)
  )

  ## Spectra numbers
  for levels in mslevelcounts.keys():
    metrics.append(
      mzqc.QualityMetric(cvRef="QC",
                         accession="QC:0000000",
                         name="Number of MS{l} spectra".format(l=str(levels)),
                         value=mslevelcounts[levels])
    )

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Number of chromatograms",
                       value=len(exp.getChromatograms()))
  )

  ## Ranges
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="MZ aquisition range",
                       value=[min_mz, max_mz])
  )

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="RT aquisition range",
                       value=[exp[0].getRT(), exp[exp.size() - 1].getRT()])
  )

  # TIC
  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Total ion current",
                       value=tic_tab)
  )

  # Chrom
  chrom_tab: Dict[str, List[Any]] = defaultdict(list)
  chroms = exp.getChromatograms()
  for t in chroms:
    if t.getChromatogramType() == oms.ChromatogramSettings.ChromatogramType.TOTAL_ION_CURRENT_CHROMATOGRAM:
      for chro_peak in t:
        chrom_tab['RT'].append(chro_peak.getRT())
        chrom_tab['int'].append(chro_peak.getIntensity())
      break

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Chromatogram",
                       value=chrom_tab)
  )
  # TODO is there a difference between TIC as defined in MS:1000235 and the chromatogram you get from TRP?? In MZML it says its a MS:1000235 (ion current detected in each of a series of mass spectra) but is it?
  # TODO consider collection of spectrum_native_id
  return mzqc.RunQuality(metadata=meta, qualityMetrics=metrics)


def describeMSCollectionTime(trap_metrics: mzqc.QualityMetric, ms_level: int) -> List[mzqc.QualityMetric]:
  """
    describeMSCollectionTime calculates the descriptive statistics metrics for ion collection times of spectra from a given level.

    From the proto-metrics on ion collection for a given MS level, the function calculates descriptive statistics metrics for
    the distribution of ion collection times from all involved mass spectra.
    Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    trap_metrics : mzqc.QualityMetric
        The proto-metrics on ion collection times from the respective MS level containing 'traptime' values.
    ms_level : int
        The MS level considered to produce the right QC metric accession

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics
    """
  metrics: List[mzqc.QualityMetric] = list()
  arr = np.array(trap_metrics['traptime'])
  q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Q1, Q2, Q3 for MS level {ms_level} trap time collection".format(ms_level=ms_level),
                       value=[q1, q2, q3])
  )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Sigma for MS level {ms_level} trap time collection".format(ms_level=ms_level),
                                    value=s)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Mean of frequency for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=m)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Frequency for MS level {ms_level} collection +/-1.5*IQR outlier".format(
                                      ms_level=ms_level),
                                    value=ol)
                 )

  return metrics


def getESIstability(ion_intensity_metric: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
  """
    getESIstability calculates the count of signal jumps and falls during the course of a mass-spectrometry run's acquisition time.

    Counts the number of signal jumps/falls of at least 10-fold magnitude.

    Parameters
    ----------
    ion_intensity_metric : mzqc.QualityMetric
        Proto-metric containing the 'int' values of signal intensity in timely order

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  folds = np.true_divide(ion_intensity_metric.value['int'][:-1], ion_intensity_metric.value['int'][1:])
  jumps = len(np.where(folds > 10)[0])
  falls = len(np.where(folds < 1 / 10)[0])

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="signal jump (10x) count",
                       value=jumps)
  )

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="signal fall (10x) count",
                       value=falls)
  )
  return metrics


def describeMSfrequency(spectrum_acquisition_metrics_MS: mzqc.QualityMetric, start_time: datetime.datetime,
                        ms_level: int) -> List[mzqc.QualityMetric]:
  """
    describeMSfrequency calculates the descriptive statistics metrics for spectra acquisition from a given level.

    From the proto-metrics on spectrum acquisition for a given MS level, the function calculates descriptive statistics metrics for
    the distribution of spectra acquisition frequency from all involved mass spectra.
    Namely, min and max, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    spectrum_acquisition_metrics_MS : mzqc.QualityMetric
        Proto-metric containing 'RT' values for all involved spectra
    start_time : datetime.datetime
        MS run start time
    ms_level : int
        The MS level considered to produce the right QC metric accession

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  rts = [start_time + datetime.timedelta(seconds=i) for i in spectrum_acquisition_metrics_MS.value['RT']]
  arr = np.array([len(list(g)) for k, g in itertools.groupby(rts, key=lambda d: d.minute)])
  q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="fastest frequency for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=max(arr))
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="slowest frequency for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=min(arr))
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Q1, Q2, Q3 of frequency for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=[q1, q2, q3])
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Sigma of frequency for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=s)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Mean of frequency for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=m)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Frequency for MS level {ms_level} collection +/-1.5*IQR outlier".format(
                                      ms_level=ms_level),
                                    value=ol)
                 )

  return metrics


def describeMSdensity(spectrum_acquisition_metrics_MS: mzqc.QualityMetric, start_time: datetime.datetime,
                      ms_level: int) -> List[mzqc.QualityMetric]:
  """
    describeMSdensity calculates the descriptive statistics metrics for spectra's peak density from a given level.

    From the proto-metrics on spectrum acquisition for a given MS level, the function calculates descriptive statistics metrics for
    the distribution of peak density from all involved mass spectra.
    Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    spectrum_acquisition_metrics_MS : mzqc.QualityMetric
        Proto-metric containing 'RT' and 'peakcount' values for all involved spectra
    start_time : datetime.datetime
        MS run start time
    ms_level : int
        The MS level considered to produce the right QC metric accession

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  rts = [start_time + datetime.timedelta(seconds=i) for i in spectrum_acquisition_metrics_MS.value['RT']]
  arr = np.array(spectrum_acquisition_metrics_MS.value['peakcount'])
  q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Q1, Q2, Q3 of peak density for MS level {ms_level} collection".format(ms_level=ms_level),
                       value=[q1, q2, q3])
  )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Sigma of peak density for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=s)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Mean of peak density for MS level {ms_level} collection".format(
                                      ms_level=ms_level),
                                    value=m)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Peak density for MS level {ms_level} collection +/-1.5*IQR outlier".format(
                                      ms_level=ms_level),
                                    value=ol)
                 )

  return metrics


def getAnalysedSignalMetrics(tandem_spectrum_metrics_MS2: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
  """
    getAnalysedSignalMetrics calculates a metric on the proportion of signal analysed with subsequent tandem spectra.

    The function calculates the median ratio of max survey scan intensity over sampled precursor intensity for the bottom (by MS1 max) half of MS2.

    Parameters
    ----------
    tandem_spectrum_metrics_MS2 : mzqc.QualityMetric
        Proto-metric of tandem spectra containing values for 'RT', 'precursor_mz', 'precursor_intensity', 'surveyscan_intensity_sum', 'surveyscan_intensity_max'.

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  # Fraction of total MS2 scans identified in the first quartile of peptides sorted by MS1 intensity (sum)
  np_prec = np.array([tandem_spectrum_metrics_MS2.value['RT'],
                      tandem_spectrum_metrics_MS2.value['precursor_mz'],
                      tandem_spectrum_metrics_MS2.value['precursor_intensity'],
                      tandem_spectrum_metrics_MS2.value['surveyscan_intensity_sum'],
                      tandem_spectrum_metrics_MS2.value['surveyscan_intensity_max']])

  # DS-3B reimpl.: median( (surv max / prec int) for bottom 50% of all precursors )
  np_prec = np_prec[:, np_prec[4].argsort()]
  # Ratio of MS1 maximum to MS1 value at sampling for bottom 50% of analytes by MS1 maximum intensity (1 = sampled at peak maxima)
  bottom_sampled_prec = np_prec[:, np_prec[4] < np.median(np_prec[4])]
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Median ratio of max survey scan intensity over sampled precursor intensity for the bottom (by MS1 max) half of MS2",
                                    value=np.median(bottom_sampled_prec[4] / bottom_sampled_prec[2]))
                 )

  return metrics


def describePrecursorIntensity(tandem_spectrum_metrics_MS2: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
  """
    describePrecursorIntensity calculates the descriptive statistics metrics for spectra's peak density from a given level.

    From the proto-metrics on tandem spectra, the function calculates descriptive statistics metrics for
    the distribution of precursor intensity. Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    tandem_spectrum_metrics_MS2 : mzqc.QualityMetric
        Proto-metric of tandem spectra containing values for 'precursor_intensity'

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  arr = np.array(tandem_spectrum_metrics_MS2.value['precursor_intensity'])
  q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Maximum precursor intensity",
                                    value=max(arr))
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Minmum precursor intensity",
                                    value=min(arr))
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Q1, Q2, Q3 of precursor intensities",
                                    value=[q1, q2, q3])
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Sigma of precursor intensities",
                                    value=s)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Mean of precursor intensities",
                                    value=m)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Precursor intensity +/-1.5*IQR outlier",
                                    value=ol)
                 )

  return metrics
