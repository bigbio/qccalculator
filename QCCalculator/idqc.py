import datetime
import itertools
import re
import sys
import warnings
import zipfile
import urllib
from collections import defaultdict
from itertools import chain
from os.path import basename
from statistics import mean, median, stdev
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

from QCCalculator import utils, noiseqc

import numpy as np
import pandas
import pronto
import pyopenms as oms
from pyopenms import ActivationMethod
from Bio import SeqIO, SeqRecord
from Bio.SeqUtils import ProtParam
from mzqc import MZQCFile as mzqc

"""
Methods to calculate general identification related quality metrics
"""

# TODO if unfiltered, filter first? If no decoy then what?
def getIDQuality(exp: oms.MSExperiment, pro_ids: List[oms.ProteinIdentification], pep_ids: List[oms.PeptideIdentification], ms2num: int = 0) -> List[mzqc.QualityMetric]:
    """
    getIDQuality calculates the id-based QualityMetrics from a mass spectrometry peak file and associated identification file. 

    Calculated are the id-based QC metrics and proto-metrics necessary to calculate more elaborate QC metrics with even more additional data (e.g. multiple runs).

    Parameters
    ----------
    exp : oms.MSExperiment
        The mass spectrometry peak file to calculate metrics from
    pro_ids : List[oms.ProteinIdentification]
        List of PyOpenMS ProteinIdentification as from reading a common identification file
    pep_ids : List[oms.PeptideIdentification]
        List of PyOpenMS PeptideIdentification as from reading a common identification file
    ms2num : int, optional
        The total number of tandem spectra as from the id-free metrics, by default 0

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """    
    metrics: List[mzqc.QualityMetric] = list() 
    params = pro_ids[0].getSearchParameters()
    # var_mods = params.variable_modifications

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sequence database name", 
                value=pro_ids[0].getSearchParameters().db)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sequence database version", 
                value=pro_ids[0].getSearchParameters().db_version)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sequence database taxonomy", 
                value=pro_ids[0].getSearchParameters().taxonomy)
    )

    spectrum_count: int = 0
    psm_count: int = 0
    runs_coun: int = 0
    protein_evidence_count: int = 0

    # TODO call mc functions
    missedcleavages: int = 0
    missedcleavages_total: int = 0

    peptides_allhits: Set[str] = set()
    peptides: Set[str] = set()
    proteins: Set[str] = set()

    for pepi in pep_ids:
      if not pepi.empty():
        # TODO if not decoy and not under threshold
        spectrum_count += 1
        psm_count += len(pepi.getHits())
        for psm in pepi.getHits():
            peptides_allhits.add(psm.getSequence().toString())
        if pepi.getHits():
            peptides.add(pepi.getHits()[0].getSequence().toString())

    for proid in pro_ids:
        protein_evidence_count += len(proid.getHits())
        for p in proid.getHits():
            proteins.add(p.getAccession())

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number of protein evidences", 
                value=protein_evidence_count)
    )

    # TODO not yet factoring in protein inference, one psm might still account for several evidences
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number of identified proteins", 
                value=len(proteins))
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number of PSM", 
                value=psm_count)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number of peptide spectra", 
                value=spectrum_count)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number identified unique peptide sequences", 
                value=len(peptides))
    )

    identification_accuracy_metrics: Dict[str,List[Any]] = defaultdict(list)
    identification_scoring_metrics: Dict[str,List[Any]] = defaultdict(list)
    identification_sequence_metrics: Dict[str,List[Any]] = defaultdict(list)
    hydrophobicity_metrics: Dict[str,List[Any]] = defaultdict(list)
    
    # TODO constants available since 2.5 as oms.Constants.PROTON_MASS_U
    # PROTON_MASS_U = 1.00727646677  # Constants::PROTON_MASS_U unavailable

    score_type = pep_ids[0].getScoreType()

    psims = utils.obtainOntology("psi-ms")

    name_indexed = {psims[x].name: psims[x] for x in psims}
    score_indexed = {x.name: x for x in chain(psims['MS:1001143'].subclasses(),psims['MS:1001153'].subclasses(),psims['MS:1002347'].subclasses(),psims['MS:1002363'].subclasses())}

    if score_type in name_indexed:
        if not score_type in score_indexed:
            warnings.warn("Score type does not correspond to a score type in the OBO, proceed at own risk.", Warning)
            score_col_name = name_indexed[score_type].id
        else:
            score_col_name = score_indexed[score_type].id
    else:
        warnings.warn("OBO does not contain any entry matching the identification score, proceed at own risk.", Warning)
        score_col_name = score_type

    for pepi in pep_ids:            
        pid = utils.pep_native_id(pepi)
        if pepi.getHits():
            tmp = pepi.getHits()[0]  # TODO apply custom filters and also consider 'pass_threshold'
            identification_scoring_metrics['RT'].append(pepi.getRT())
            identification_scoring_metrics['c'].append(tmp.getCharge())
            identification_scoring_metrics[score_col_name].append(tmp.getScore())

            tw = (tmp.getSequence().getMonoWeight(0,0) + tmp.getCharge() * oms.Constants.PROTON_MASS_U) / tmp.getCharge()
            dppm = utils.getMassDifference(tw, pepi.getMZ(), True)      
            identification_accuracy_metrics['RT'].append(pepi.getRT())
            identification_accuracy_metrics['MZ'].append(pepi.getMZ())
            identification_accuracy_metrics['delta_ppm'].append(dppm)
            err = utils.getMassDifference(tw, pepi.getMZ(), False) 
            identification_accuracy_metrics['abs_error'].append(err)

            hydrophobicity_metrics['RT'].append(pepi.getRT())
            hydrophobicity_metrics['gravy'].append(ProtParam.ProteinAnalysis(tmp.getSequence().toUnmodifiedString()).gravy())

            identification_sequence_metrics['RT'].append(pepi.getRT())
            identification_sequence_metrics['peptide'].append(tmp.getSequence().toString().lstrip().rstrip())
            identification_sequence_metrics['target'].append(tmp.getMetaValue('target_decoy').lower() == 'target')
            identification_sequence_metrics['native_id'].append(pid)

    #   #varmod???         
    #   for (UInt w = 0; w < var_mods.size(); ++w)
    #   {
    #     at.colTypes.push_back(String(var_mods[w]).substitute(' ', '_'));
    #   }

    ## Basic id metrics
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identification scoring metric values", 
                value=identification_scoring_metrics)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identifications accuracy metric values", 
                value=identification_accuracy_metrics)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Hydrophobicity metric values", 
                value=hydrophobicity_metrics)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identifications sequence metric values", 
                value=identification_sequence_metrics)
    )

    ## simple id metrics 
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identification to tandem spectra ratio", 
                value=float(len(pep_ids))/float(ms2num))
    )

    return metrics

def describeChargeRatios(identification_scoring_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    describeChargeRatios calculates the descriptive statistics metrics for charge ratios of identified tandem spectra.

    From the proto-metrics on identification scores, the function calculates descriptive statistics metrics on the 
    charge ratios from all identified tandem spectra. Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    identification_scoring_metrics : mzqc.QualityMetric
        The proto-metrics on identification scores containing 'c' values.

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics 
    """    
    metrics: List[mzqc.QualityMetric] = list() 
    if 'c' not in identification_scoring_metrics:
        warnings.warn("No charges in given annotation, ignoring charge ratio metrics.", Warning)
        return metrics

    # IS3A <- c1n / c2n
    # IS3B <- c3n / c2n
    # IS3C <- c4n / c2n
    unique_charges, charge_freq = np.unique(identification_scoring_metrics.value['c'], return_counts=True)
    c1i = np.where(unique_charges == 1)
    c1n = charge_freq[c1i[0][0]] if len(c1i[0])>0 else 0 
    c2i = np.where(unique_charges == 2)
    c2n = charge_freq[c2i[0][0]] if len(c2i[0])>0 else 0 
    c3i = np.where(unique_charges == 3)
    c3n = charge_freq[c3i[0][0]] if len(c3i[0])>0 else 0 
    c4i = np.where(unique_charges == 4)
    c4n = charge_freq[c4i[0][0]] if len(c4i[0])>0 else 0 
    mi = c4i if c4i > 0 else c3i if c3i > 0 else c2i if c2i > 0 else c1i if c1i > 0 else 0
    c5p = sum(charge_freq[mi+1:])

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="IS3A", 
                value= (c1n / c2n) if c2n > 0 else np.NAN )
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="IS3B", 
                value= (c3n / c2n) if c2n > 0 else np.NAN )
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="IS3C", 
                value= (c4n / c2n) if c2n > 0 else np.NAN )
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="IS3X", 
                value=np.median(identification_scoring_metrics.value['c']))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="IS3Y", 
                value=np.mean(identification_scoring_metrics.value['c']))
    )

    return metrics

def describeErrorRates(identification_accuracy_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    describeErrorRates calculates the descriptive statistics metrics for charge ratios of identified tandem spectra.

    From the proto-metrics on identification accuracy, the function calculates descriptive statistics metrics on the 
    error rates from all identified tandem spectra. Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    identification_accuracy_metrics : mzqc.QualityMetric
        The proto-metrics on identification accuracies containing 'delta_ppm' and 'abs_error' values.

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics 
    """    
    metrics: List[mzqc.QualityMetric] = list() 
    if 'delta_ppm' not in identification_accuracy_metrics:
        warnings.warn("No error values in given annotation, ignoring identification error rate metrics.", Warning)
        return metrics

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15A", 
                value= np.median(identification_accuracy_metrics.value['abs_error']) )
    )
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15B", 
                value=np.mean(identification_accuracy_metrics.value['abs_error']) )
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15C", 
                value=np.median(identification_accuracy_metrics.value['delta_ppm']) )
    )

    arr = np.array(identification_accuracy_metrics.value['delta_ppm'])
    q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15D", 
                value=q3-q1)
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm Q1, Q2, Q3", 
                value=[q1,q2,q3])
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm sigma", 
                value=s)
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm mean", 
                value=m)
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm +/-1.5*IQR outlier", 
                value=ol)
    )

    return metrics

# TODO target decoy separation and FDR, correct name would be getSamplingRates not ratios
def getSamplingRatios(identification_sequence_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    getSamplingRatios calculates the sampling ratio metric for identified tandem spectra.

    From the proto-metrics on identified sequences, the function calculates sampling rate and frequency.

    Parameters
    ----------
    identification_sequence_metrics : mzqc.QualityMetric
        The proto-metrics on identified sequences containing 'peptide' (sequence) values.

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics 
    """    
    metrics: List[mzqc.QualityMetric] = list() 

    sample_rate, sample_rate_counts = np.unique(np.unique(identification_sequence_metrics.value['peptide'], return_counts=True)[1], return_counts=True)
    # explicitly enum all sampling rates up to the max.
    explicit_rate_counts = np.zeros( np.max(sample_rate) ) 
    explicit_rate = np.arange(1, np.max(sample_rate)+1)
    explicit_indices = np.where(np.isin(explicit_rate,sample_rate))
    np.put(explicit_rate,explicit_indices,sample_rate)
    np.put(explicit_rate_counts,explicit_indices,sample_rate_counts)
    
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sampling frequencies", 
                value={'sampling rate': list(explicit_rate), 
                    'frequencies': list(explicit_rate_counts)})
    )
      
    return metrics

def describeIdentificationScores(identification_scoring_metrics:mzqc.QualityMetric, score_type:str) -> List[mzqc.QualityMetric]:
    """
    describeIdentificationScores calculates the descriptive statistics metrics for the scoring of identified tandem spectra.

    From the proto-metrics on identification scores, the function calculates descriptive statistics metrics on the 
    charge id scores from all identified tandem spectra. Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    identification_scoring_metrics : mzqc.QualityMetric
        The proto-metrics on identification scores containing `score_type` values.
    score_type : str
        The score_type descriptor used to create the identification score values category from `identification_scoring_metrics`

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics 
    """    
    metrics: List[mzqc.QualityMetric] = list()    
    qs =  np.quantile(identification_scoring_metrics.value[score_type], [.25,.5,.75])

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="scores Q1, Q2, Q3", 
                value=list(qs))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="score sigma", 
                value=np.std(identification_scoring_metrics.value[score_type]))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="score mean", 
                value=np.mean(identification_scoring_metrics.value[score_type]))
    )

    np_score = np.array(identification_scoring_metrics.value[score_type])
    low_out = qs[0]-(1.5*(qs[2]-qs[0]))
    high_out = qs[2]+(1.5*(qs[2]-qs[0]))
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="score +/-1.5*IQR outlier", 
                value=np.extract((np_score<low_out) | (np_score>high_out), np_score))
    )

    return metrics

def getIdentifiedSignalMetrics(tandem_spectrum_metrics_MS2:mzqc.QualityMetric, 
        spectrum_acquisition_metrics_MS1: mzqc.QualityMetric,
        identification_accuracy_metrics: mzqc.QualityMetric,
        tic_table: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    getIdentifiedSignalMetrics calculate metrics on the proportions of recorded signal identified.

    The metrics calculated include the median ratio of max survey scan intensity over sampled precursor 
    intensity for peptides identified, the fractions of identified MS2 in precursor intensity Quartiles, 
    median SN for MS1 spectra in RT range in which the first half of peptides are identified, and 
    median TIC value of RT range in which half of peptides are identified.

    Parameters
    ----------
    tandem_spectrum_metrics_MS2 : mzqc.QualityMetric
        The proto-metrics on tandem spectra containing 'RT', 'precursor_mz', 'precursor_intensity', 'surveyscan_intensity_sum', 'surveyscan_intensity_max' values.
    spectrum_acquisition_metrics_MS1 : mzqc.QualityMetric
        The proto-metrics on MS1 spectra containing 'RT' and 'SN' values
    identification_accuracy_metrics : mzqc.QualityMetric
        The proto-metrics on identification accuracies containing 'RT' and 'MZ' values
    tic_table : mzqc.QualityMetric
        The proto-metrics on total ion current intensities containing 'RT' and 'int'

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics 
    """    
    metrics: List[mzqc.QualityMetric] = list() 

    # Fraction of total MS2 scans identified in the first quartile of peptides sorted by MS1 intensity (sum)
    np_prec = np.array([tandem_spectrum_metrics_MS2.value['RT'], 
                        tandem_spectrum_metrics_MS2.value['precursor_mz'], 
                        tandem_spectrum_metrics_MS2.value['precursor_intensity'], 
                        tandem_spectrum_metrics_MS2.value['surveyscan_intensity_sum'], 
                        tandem_spectrum_metrics_MS2.value['surveyscan_intensity_max']])

    # DS-3A reimpl.: median( (surv max / prec int) for all ident. prec ) 
    id_coord = np.array([identification_accuracy_metrics.value['RT'],identification_accuracy_metrics.value['MZ']])  # TODO make sure intersection is round-proof
    intersected = np.intersect1d(np_prec[1],id_coord[1], assume_unique=False, return_indices=True)
    np_id = np_prec[:,intersected[1]]
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median ratio of max survey scan intensity over sampled precursor intensity for peptides identified", 
                value=np.median(np_id[4] / np_id[2]))
    )

    # MS1-3A reimpl.: Ratio of 95th over 5th percentile MS1 maximum intensity values for identified peptides (approximates dynamic range of signal)
    p05, p95 = np.quantile(np_id[4], [.05, .95])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Ratio of 95th over 5th percentile of precursor intensity for identified peptides", 
                value=p95 / p05 )
    )

    # Quartiles by MS1 maximum intensity
    # Fraction of identified MS2 Spectra within
    # MS2-4A : 0 and Q1
    # MS2-4B : Q1 and Q2
    # MS2-4C : Q2 and Q3
    # MS2-4D : above Q3
    q1,q2,q3 = np.quantile(np_prec[4], [.25, .5, .75])
    
    tandem_upto_q1 = np.shape(np_prec[:,np_prec[4]<q1])[1]
    id_upto_q1 = np.shape(np_id[:,np_id[4]<q1])[1]

    tandem_between_q1q2 = np.shape(np_prec[:,(q1<np_prec[4]) & (np_prec[4]<q2)])[1]
    id_between_q1q2 = np.shape(np_id[:,(q1<np_id[4]) & (np_id[4]<q2)])[1]

    tandem_between_q2q3 = np.shape(np_prec[:,(q2<np_prec[4]) & (np_prec[4]<q3)])[1]
    id_between_q2q3 = np.shape(np_id[:,(q2<np_id[4]) & (np_id[4]<q3)])[1]

    tandem_above_q3 = np.shape(np_prec[:,q3<np_prec[4]])[1]
    id_above_q3 = np.shape(np_id[:,q3<np_id[4]])[1]
    
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Fraction of identified MS2 below Q1 of precursor intensity.", 
                value=tandem_upto_q1 / id_upto_q1)
    )
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Fraction of identified MS2 between Q1 and Q2 of precursor intensity.", 
                value=tandem_between_q1q2 / id_between_q1q2)
    )
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Fraction of identified MS2 between Q2 and Q3 of precursor intensity.", 
                value=tandem_between_q2q3 / id_between_q2q3)
    )
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Fraction of identified MS2 above Q3 of precursor intensity.", 
                value=tandem_above_q3 / id_above_q3)
    )

    # MS1-3B reimpl.: Median maximum MS1 value for identified peptides
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median precursor intensity of identified MS2", 
                value=np.median(np_id[4]))
    )

    # MS1-2A Median SN for MS1 spectra in RT range in which half (which half???) of peptides are identified	
    np_id = np_id[:,np_id[0].argsort()]
    msn = np.array([spectrum_acquisition_metrics_MS1.value['RT'], spectrum_acquisition_metrics_MS1.value['SN']])	
    median_id_rt = np.quantile(np_id[0], [.5])[0]
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median SN for MS1 spectra in RT range in which the first half of peptides are identified", 
                value=np.median(msn[:, msn[0]<median_id_rt ][1]) )
    )

    # the smallest rt range which contains half of all identified Spectra
    half_id_size = np.round(np_id.shape[1]/2)
    all_diff = np_id[:,-1*int(half_id_size):][0] - np_id[:,:int(half_id_size)][0]  # the last (half_id_size) many - the first (half_id_size) many
    min_start_index = np.argmin(all_diff)
    min_stop_index = min_start_index+int(half_id_size)-1
    rt_interval = np_id[0,min_stop_index] - np_id[0,min_start_index]
        
    densest = np.median(msn[:,(np_id[0,min_start_index]<msn[0]) & (msn[0]<np_id[0,min_stop_index])][1])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median SN for MS1 spectra in densest RT range in which any half of peptides are identified", 
                value=densest)
    )

    # is median tic value really meaningful? My guess is ratio of tic sum of RT half identified and rest is a better indicator (>1: most signal is in the most exlained region)
    # Median TIC value of RT range in which half of peptides are identified	
    np_tic = np.array([tic_table.value['RT'], tic_table.value['int']])
    densest_id_tic = np.median(np_tic[:,(np_id[0,min_start_index]<np_tic[0]) & (np_tic[0]<np_id[0,min_stop_index])][1])
    
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median TIC value of RT range in which half of peptides are identified", 
                value=np.median(np_tic[:, np_tic[0]<median_id_rt ][1]) )
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median TIC value of densest RT range in which any half of peptides are identified", 
                value=densest_id_tic )
    )

    return metrics

def describeIdentifiedPrecursorIntensity(tandem_spectrum_metrics_MS2:mzqc.QualityMetric, 
        identification_accuracy_metrics: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    """
    describeIdentifiedPrecursorIntensity calculates the descriptive statistics metrics for precursor intensities of identified tandem spectra.

    From the proto-metrics on identification accuracies and tandem spectra, the function calculates descriptive statistics metrics on the 
    precursor intensities from all identified tandem spectra. Namely, min and max, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    tandem_spectrum_metrics_MS2 : mzqc.QualityMetric
        The proto-metrics on tandem spectra containing 'RT', 'precursor_mz', 'precursor_intensity', 'surveyscan_intensity_sum', 'surveyscan_intensity_max' values.
    identification_accuracy_metrics : mzqc.QualityMetric
        The proto-metrics on identification accuracies containing 'RT' and 'MZ' values

    Returns
    -------
    List[mzqc.QualityMetric]
        The list of metrics 
    """    
    metrics: List[mzqc.QualityMetric] = list() 

    # Fraction of total MS2 scans identified in the first quartile of peptides sorted by MS1 intensity (sum)
    np_prec = np.array([tandem_spectrum_metrics_MS2.value['RT'], 
                        tandem_spectrum_metrics_MS2.value['precursor_mz'], 
                        tandem_spectrum_metrics_MS2.value['precursor_intensity'], 
                        tandem_spectrum_metrics_MS2.value['surveyscan_intensity_sum'], 
                        tandem_spectrum_metrics_MS2.value['surveyscan_intensity_max']])

    # DS-3A reimpl.: median( (surv max / prec int) for all ident. prec ) 
    id_coord = np.array([identification_accuracy_metrics.value['RT'],identification_accuracy_metrics.value['MZ']])  # TODO make sure intersection is round-proof
    intersected = np.intersect1d(np_prec[1],id_coord[1], assume_unique=False, return_indices=True)
    np_id = np_prec[:,intersected[1]]
    arr = np_id[2]
    q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Maximum identified precursor intensity", 
                value=max(arr))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Minmum identified precursor intensity", 
                value=min(arr))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Q1, Q2, Q3 of identified precursor intensities", 
                value=[q1, q2, q3])
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sigma of identified precursor intensities",
                value=s)
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Mean of identified precursor intensities",
                value=m)
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Precursor identified intensity +/-1.5*IQR outlier",
                value=ol)
    )

    return metrics
