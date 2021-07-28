import datetime
import itertools
import numpy as np
import pandas as pd
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import pyopenms as oms
import configparser
import scipy.ndimage as nd

from mzqc import MZQCFile as mzqc
from qccalculator import utils

def calcESIStability(ion_intensity_metric:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
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

    folds = np.true_divide(ion_intensity_metric.value['int'][:-1],ion_intensity_metric.value['int'][1:])
    jumps = len(np.where(folds > 10)[0])
    falls = len(np.where(folds < 1/10)[0])

    metrics.append(
        mzqc.QualityMetric(accession="QC:4000142", 
                name="signal jump (10x) count", 
                value=jumps)
    )

    metrics.append(
        mzqc.QualityMetric(accession="QC:4000143", 
                name="signal fall (10x) count", 
                value=falls)
    )
    return metrics

def calcMSFrequency(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, start_time: datetime.datetime, ms_level: int) -> List[mzqc.QualityMetric]:
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
    arr = np.array([len(list(g)) for k, g in itertools.groupby(rts, key=lambda d: d.minute)])  # number of spectra per each minute of acquisition

    quarter_arrs = np.array_split(arr,4)
    hz_per_quarter = [q.sum()/(q.size*60) for q in quarter_arrs]

    if ms_level == 1:
        acfmax = "QC:4000115"
        acfmin = "QC:4000117"
        acrtq = "QC:4000120"
    elif ms_level == 2:
        acfmax = "QC:4000116"
        acfmin = "QC:4000118"
        acrtq = "QC:4000121"
    else:
        acfmax = "QC:0000000"
        acfmin = "QC:0000000"
        acrtq = "QC:0000000"

    metrics.append(mzqc.QualityMetric(accession=acfmax, 
                name="Fastest frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=max(arr))
    )

    metrics.append(mzqc.QualityMetric(accession=acfmin, 
                name="Slowest frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=min(arr))
    )

    metrics.append(mzqc.QualityMetric(accession=acfmin, 
                name="MS{ms_level} frequency in RT quarters".format(ms_level=ms_level), 
                value=hz_per_quarter)
    )


    return metrics

def calcMSDensity(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, start_time: datetime.datetime, ms_level: int) -> List[mzqc.QualityMetric]:
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

    if ms_level == 1:
        acm = "QC:4000131"
        acs = "QC:4000132"
        acq = "QC:4000130"
        aco = "QC:4000133"
    elif ms_level == 2:
        acm = "QC:4000136"
        acs = "QC:4000136"
        acq = "QC:4000134"
        aco = "QC:4000137"
    else:
        acm = "QC:0000000"
        acs = "QC:0000000"
        acq = "QC:0000000"
        aco = "QC:0000000"

    metrics.append(
        mzqc.QualityMetric(accession=acq, 
                name="Peak density distribution MS{ms_level} Q1, Q2, Q3".format(ms_level=ms_level), 
                value=[q1, q2, q3])
    )

    metrics.append(mzqc.QualityMetric(accession=acs, 
                name="Peak density distribution MS{ms_level} sigma".format(ms_level=ms_level),
                value=s)
    )

    metrics.append(mzqc.QualityMetric(accession=acm, 
                name="Peak density distribution MS{ms_level} mean".format(ms_level=ms_level),
                value=m)
    )

    metrics.append(mzqc.QualityMetric(accession=aco, 
                name="Peak density distribution MS{ms_level} outliers (<Q1-1.5*IQR, >Q3+1.5*IQR)".format(ms_level=ms_level),
                value=ol)
    )

    return metrics

def calcAnalysedSignalMetrics(tandem_spectrum_metrics:mzqc.QualityMetric, isolation_window_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    getAnalysedSignalMetrics calculates a metric on the proportion of signal analysed with subsequent tandem spectra.

    The function calculates the median ratio of max survey scan intensity over sampled precursor intensity for the bottom (by MS1 max) half of MS2.

    Parameters
    ----------
    tandem_spectrum_metrics : mzqc.QualityMetric
        Proto-metric of tandem spectra containing values for RT, precursor_mz, precursor_intensity, surveyscan_intensity_sum, surveyscan_intensity_max.
    isolation_window_metrics : mzqc.QualityMetric
        Proto-metric of isolation windows and precursor intensities containing values for RT, precursor_mz, precursor_intensity, surveyscan_intensity_sum, surveyscan_intensity_max.

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """    
    metrics: List[mzqc.QualityMetric] = list() 

    np_prec =  pd.merge(pd.DataFrame(tandem_spectrum_metrics.value),
                        pd.DataFrame(isolation_window_metrics.value), 
                        how="inner", on='RT')

    # Fraction of total MS2 scans identified in the first quartile of peptides sorted by MS1 intensity (sum)

    # DS-3B reimpl.: median( (surv max / prec int) for bottom 50% of all precursors ) 
    # Ratio of MS1 maximum to MS1 value at sampling for bottom 50% of analytes (higher is sampled farther from the apex)
    np_prec.sort_values(by=['precursor_int'])
    bottom_sampled_prec = np_prec[np_prec.precursor_int < np_prec.precursor_int.median()]
    
    metrics.append(mzqc.QualityMetric(accession="QC:4000108", 
                name="Explained base peak intensity median", 
                value=np.median(bottom_sampled_prec.surveyscan_int_max / bottom_sampled_prec.precursor_int))
    )

    return metrics

def calcPrecursorIntensityMetrics(tandem_spectrum_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    describePrecursorIntensity calculates the descriptive statistics metrics for spectra's peak density from a given level.

    From the proto-metrics on tandem spectra, the function calculates descriptive statistics metrics for 
    the distribution of precursor intensity. Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    tandem_spectrum_metrics : mzqc.QualityMetric
        Proto-metric of tandem spectra containing values for 'precursor_int'

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """    
    metrics: List[mzqc.QualityMetric] = list() 
    
    arr = np.array(tandem_spectrum_metrics.value['precursor_int'])
    q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

    metrics.append(mzqc.QualityMetric(accession="QC:4000119", 
                name="Precursor intensity range", 
                value=[min(arr),max(arr)])
    )

    metrics.append(mzqc.QualityMetric(accession="QC:4000138", 
                name="Precursor intensity distribution Q1, Q2, Q3", 
                value=[q1, q2, q3])
    )

    metrics.append(mzqc.QualityMetric(accession="QC:4000140", 
                name="Precursor intensity distribution sigma",
                value=s)
    )

    metrics.append(mzqc.QualityMetric(accession="QC:4000139", 
                name="Precursor intensity distribution mean",
                value=m)
    )

    metrics.append(mzqc.QualityMetric(accession="QC:4000141", 
                name="Precursor intensity distribution outliers (<Q1-1.5*IQR, >Q3+1.5*IQR)",
                value=ol)
    )

    return metrics

def calcSNMetrics(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, ms_level: int) -> List[mzqc.QualityMetric]:    
    """
    calcSNMetrics collect S/N related QC metrics from a super metric collected in a first pass of the input mzML

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
    q1, q2, q3, s, m, ol = utils.extractDistributionStats(np.array(spectrum_acquisition_metrics_MS.value['SN']))

    if ms_level == 1: 
        acq = "QC:4000184"
        acm = "QC:4000185"
        aco = "QC:4000187"
        acs = "QC:4000186"
    elif ms_level == 2: 
        acq = "QC:4000188"
        acm = "QC:4000189"
        aco = "QC:4000191"
        acs = "QC:4000190"
    else:
        acq = "QC:0000000"
        acm = "QC:0000000"
        aco = "QC:0000000"
        acs = "QC:0000000"

    metrics.append(mzqc.QualityMetric(
                accession=acq,
                name="Signal-to-noise ratio in MS{ms_level} - Q1, Q2, Q3".format(ms_level=ms_level),
                value=list([q1,q2,q3]))
    )

    metrics.append(mzqc.QualityMetric( 
                accession=acs, 
                name="Signal-to-noise ratio in MS{ms_level} - sigma".format(ms_level=ms_level),
                value=s)
    )

    metrics.append(mzqc.QualityMetric( 
                accession=acm, 
                name="Signal-to-noise ratio in MS{ms_level} - mean".format(ms_level=ms_level), 
                value=m)
    )

    metrics.append(mzqc.QualityMetric(
                accession=aco,  
                name="Signal-to-noise ratio in MS{ms_level} - outliers (<Q1-1.5*IQR, >Q3+1.5*IQR)".format(ms_level=ms_level),
                value=ol)
    )
      
    return metrics

def calcTrapTimeMetrics(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, ms_level: int) -> List[mzqc.QualityMetric]:    
    """
    calcTrapTimeMetrics collect (ion) trap time related QC metrics from a super metric collected in a first pass of the input mzML

    The trap time from each spectrum are computed into 'spectrum acquisition metrics' for each MS level, from there the trap time 
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
    q1, q2, q3, s, m, ol = utils.extractDistributionStats(np.array(spectrum_acquisition_metrics_MS.value['traptime']))

    if ms_level == 1: 
        acq = "QC:4000122"
        acm = "QC:4000123"
        aco = "QC:4000125"
        acs = "QC:4000124"
    elif ms_level == 2: 
        acq = "QC:4000126"
        acm = "QC:4000127"
        aco = "QC:4000129"
        acs = "QC:4000128"
    else:
        acq = "QC:0000000"
        acm = "QC:0000000"
        aco = "QC:0000000"
        acs = "QC:0000000"

    metrics.append(mzqc.QualityMetric( 
                accession=acq,
                name="MS{ms_level} trap time Q1, Q2, Q3".format(ms_level=ms_level),
                value=list([q1,q2,q3]))
    )

    metrics.append(mzqc.QualityMetric( 
                accession=acs, 
                name="MS{ms_level} trap time sigma".format(ms_level=ms_level),
                value=s)
    )

    metrics.append(mzqc.QualityMetric( 
                accession=acm, 
                name="MS{ms_level} trap time mean".format(ms_level=ms_level), 
                value=m)
    )

    metrics.append(mzqc.QualityMetric( 
                accession=aco,  
                name="MS{ms_level} trap time outliers (<Q1-1.5*IQR, >Q3+1.5*IQR)".format(ms_level=ms_level),
                value=ol)
    )
      
    return metrics

def calcPrecursorChargeMetrics(tandem_spectrum_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    hist = np.histogram(tandem_spectrum_metrics.value['precursor_c'],bins=range(0,max(tandem_spectrum_metrics.value['precursor_c'])))
    ms2count = len(tandem_spectrum_metrics.value['precursor_c'])

    metrics.append(mzqc.QualityMetric(accession="QC:4000064", 
                name="MS2 unknown and likely precursor charges fractions",
                value=hist[0]/ms2count)
    )

    try:
        c12 = hist[0][1]/hist[0][2]
    except:
        c12 = np.nan
    metrics.append(mzqc.QualityMetric(accession="QC:4000149", 
                name="Charged spectra ratio +1 over +2",
                value=c12)
    )

    try:
        c32 = hist[0][3]/hist[0][2]
    except:
        c32 = np.nan
    metrics.append(mzqc.QualityMetric(accession="QC:4000150", 
                name="Charged spectra ratio +3 over +2",
                value=c32)
    )

    try:
        c42 = hist[0][4]/hist[0][2]
    except:
        c42 = np.nan
    metrics.append(mzqc.QualityMetric(accession="QC:4000151", 
                name="Charged spectra ratio +4 over +2",
                value=c42)
    )

    metrics.append(mzqc.QualityMetric(accession="QC:4000152", 
                name="Mean precursor charge in all MS2",
                value=np.mean(hist[0]))
    )

    metrics.append(mzqc.QualityMetric(accession="QC:4000153", 
                name="Median precursor charge in all MS2",
                value=np.median(hist[0]))
    )

    return metrics

def calcTICSpread(tic_tab: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    
    stic = np.array(tic_tab.value['int']).sum()  # quarter the tic by accumulated intensity
    npc = np.array([tic_tab.value['RT'],tic_tab.value['int'],np.cumsum(tic_tab.value['int'])])

    q1filter = npc[0:2,(npc[2,:] < stic/4)]  # 1st quarter in respect to int
    q2filter = npc[0:2,(npc[2,:] >= stic/4) & (npc[2,:] < 2*(stic/4))]  # 2nd quarter in respect to int
    q3filter = npc[0:2,(npc[2,:] >= 2*(stic/4)) & (npc[2,:] < 3*(stic/4))]  # 3rd quarter in respect to int
    q4filter = npc[0:2,(npc[2,:] >= 3*(stic/4))]  # 4th quarter in respect to int

    # last RT of each quarter (except last)
    metrics.append(mzqc.QualityMetric(accession="QC:4000054", 
                name="RT over TIC quantile",
                value=[max(q1filter[0]), max(q2filter[0]), max(q3filter[0])])
    )

    qtic = max(tic_tab.value['RT'])/4  # quarter the tic by RT
    r1filter = npc[0:2,(npc[0,:] < qtic/4)]  # 1st quarter in respect to RT
    r2filter = npc[0:2,(npc[0,:] >= qtic/4) & (npc[0,:] < 2*(qtic/4))]  # 2nd quarter in respect to RT
    r3filter = npc[0:2,(npc[0,:] >= 2*(qtic/4)) & (npc[0,:] < 3*(qtic/4))]  # 3rd quarter in respect to RT
    r4filter = npc[0:2,(npc[0,:] >= 3*(qtic/4))]  # 4th quarter in respect to RT

    # log ratio of quarter2 by quarter1, quarter3 by quarter1, quarter4 by quarter1
    metrics.append(mzqc.QualityMetric(accession="QC:4000058", 
                name="MS1 quantile TIC ratio to Q1",
                value=[np.log(r2filter[1].sum())-np.log(r1filter[1].sum()), 
                        np.log(r3filter[1].sum())-np.log(r1filter[1].sum()), 
                        np.log(r4filter[1].sum())-np.log(r1filter[1].sum())])
    )

    # TODO implement QC:4000057 ; BUT I have not an inckling of an idea what the TIC change in that context is supposed to mean
    return metrics

def calcMSSpread(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, ms_level: int) -> List[mzqc.QualityMetric]: 
    metrics: List[mzqc.QualityMetric] = list() 

    if ms_level == 1:
        ac = "QC:4000055"
    elif ms_level == 2:
        ac = "QC:4000056"

    # def: "The interval used for acquisition of the first, second, third, and fourth quarter of all MS2 events divided by RT-Duration." 
    event_quarter1, event_quarter2, event_quarter3, event_quarter4 = np.array_split(spectrum_acquisition_metrics_MS.value['RT'], 4)
    metrics.append(mzqc.QualityMetric(accession=ac, 
                name="MS{ms_level} quantiles RT fraction".format(ms_level=ms_level),
                value=[event_quarter1.size/(np.max(event_quarter1)-np.min(event_quarter1)),
                        event_quarter2.size/(np.max(event_quarter2)-np.min(event_quarter2)),
                        event_quarter3.size/(np.max(event_quarter3)-np.min(event_quarter3)),
                        event_quarter4.size/(np.max(event_quarter4)-np.min(event_quarter4)),
                    ])
    )

    return metrics

def calcMSMSPeakDominance(exp: oms.MSExperiment, config: configparser.ConfigParser, verbose: bool=False) -> mzqc.RunQuality:
    metrics: List[mzqc.QualityMetric] = list() 

    spec_q: Dict[str,list] = {'native_id': list(), 'peak_dominance': list() }

    for spec in exp:
        if spec.getMSLevel() > 1:
            spec.sortByIntensity(1)  # desc.
            spectic = np.sum(spec.get_peaks()[1])  # spectrum total intensity
            speccs = np.cumsum(spec.get_peaks()[1])  # max peaks first cumulative intensity
            spec_q['peak_dominance'].append(speccs[np.where(speccs < spectic/2)].size)  # num of peaks to explain 50%
            spec_q['native_id'].append(utils.getSpectrumNativeID(spec))

    metrics.append(mzqc.QualityMetric(accession="QC:4000268", 
                name="Min. peaks for 50% of total spectra intensity",
                value=spec_q)
    )
    return metrics

def calcBaseDrift(tic_tab: mzqc.QualityMetric) -> mzqc.RunQuality:
    metrics: List[mzqc.QualityMetric] = list() 

    # tophat for baseline drift estimate
    # savgol filter for smoothing before tophat?

    chrom = pd.DataFrame(tic_tab.value)
    blpe_factor = 0.1  # base line points estimate factors
    struct_pts = int(round(chrom.int.size*blpe_factor))

    chrom["tophat"] = nd.white_tophat(chrom.int, None, np.repeat([1], struct_pts))
    # di = np.sum(chrom.int-chrom.tophat)
    #np.sum(chrom.int)
    metrics.append(mzqc.QualityMetric(accession="QC:4000270", 
                name="Estimated baseline drift fraction of total intensity",
                value=np.sum(chrom.int-chrom.tophat)/np.sum(chrom.int))
    )

