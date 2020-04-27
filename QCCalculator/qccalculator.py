import datetime
import hashlib
import io
import itertools
import re
import sys
import urllib
import urllib.request
import warnings
import zipfile
from collections import defaultdict
from itertools import chain
from os.path import basename
from statistics import mean, median, stdev
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas
import pronto
import pyopenms as oms
from pyopenms import ActivationMethod
from Bio import SeqIO, SeqRecord
from Bio.SeqUtils import ProtParam
from mzqc import MZQCFile as mzqc


def sha256fromfile(filename: str) -> str:
    sha  = hashlib.sha256()
    b  = bytearray(128*1024)
    mv = memoryview(b)
    with open(filename, 'rb', buffering=0) as f:
        for n in iter(lambda : f.readinto(mv), 0):
            sha.update(mv[:n])
    return sha.hexdigest()

def cast_if_int(pot_int):
    try:
        return int(pot_int)
    except ValueError as e:
        return pot_int

def pep_native_id(p: oms.Peptide) -> Union[int,None]:
    spre = p.getMetaValue('spectrum_reference')
    if spre:
        matches = re.findall("scan=(\d+)$", spre.decode())
        if len(matches)!=1:  # should really never be >1 with the `$`
            return None
        else:
            return cast_if_int(matches[0])
    else:
        return None

def spec_native_id(s: oms.MSSpectrum) -> Union[int,None]:
    spre = s.getNativeID()
    if spre:
        matches = re.findall("scan=(\d+)$", spre.decode())
        if len(matches)!=1:  # should really never be >1 with the `$`
            return None
        else:
            return cast_if_int(matches[0])
    else:
        return None

def getMassDifference(theo_mz: float, exp_mz: float, use_ppm: bool=True)-> float:
    error: float = (exp_mz - theo_mz)
    if use_ppm:
        error = error / (theo_mz * 1e-6)
    return error
  
def getSN_medianmethod(spec: oms.MSSpectrum, norm: bool=True) -> float:
    if spec.size() == 0: 
        return 0.0

    median: float = 0.0
    maxi: float = 0.0
    spec.sortByIntensity(False)
    
    mar = np.array([s.getIntensity() for s in spec])
    median = np.median(mar)

    if (not norm):
        return np.max(mar) / median

    sig = np.sum(mar[mar<=median])/mar[mar<=median].size
    noi = np.sum(mar[mar>median])/mar[mar>median].size
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
          
    return sig/noi

def getTrapTime(spec: oms.MSSpectrum, acqusition_unavailable= True) -> float:
    tt = -1.0
    if acqusition_unavailable:
        if spec.metaValueExists('MS:1000927'):
            tt = spec.getMetaValue('MS:1000927')        
    else:
        # TODO incomplete AcqusitionInfo.pxd
        if not spec.getAcquisitionInfo().empty():
            for j in spec.getAcquisitionInfo():
                if j.metaValueExists("MS:1000927"):
                    tt = j.getMetaValue("MS:1000927")
                    break
    return tt

from enum import Enum, unique

def getBasicQuality(exp: oms.MSExperiment, verbose: bool=False) -> mzqc.RunQuality:
    metrics: List[mzqc.QualityMetric] = list()
    if exp.getExperimentalSettings().getSourceFiles():
        parent_base_name: str = basename(exp.getExperimentalSettings().getSourceFiles()[0].getNameOfFile().decode())
        parent_chksm: str = exp.getExperimentalSettings().getSourceFiles()[0].getChecksum().decode()
        parent_chksm_type: str = exp.getExperimentalSettings().getSourceFiles()[0].getChecksumType()
    
    input_loc: str = exp.getExperimentalSettings().getLoadedFilePath().decode()
    base_name: str = basename(input_loc)
    chksm: str = sha256fromfile(exp.getExperimentalSettings().getLoadedFilePath().decode())
    cmpltn: str = exp.getDateTime().get().decode()
    # strt:datetime.datetime = datetime.datetime.strptime(cmpltn, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(seconds=exp.getChromatograms()[0][exp.getChromatograms()[0].size()-1].getRT()*60)

    meta: mzqc.MetaDataParameters = mzqc.MetaDataParameters(
        inputFiles=[
            mzqc.InputFile(name=base_name,location=input_loc, 
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
                                value=exp.getInstrument().getName().decode()
                            )
                            # TODO integrate parent location and checksum
                            # id: MS:1002846
                            # name: Associated raw file URI
                            # fitting checksum cv missing
                        ]
            )
        ], 
        analysisSoftware=[
            mzqc.AnalysisSoftware(cvRef="MS", accession="MS:1000752", name="TOPP software", version="2.4", uri="openms.de")
        ]
    )

    # this is mighty important to sort by RT
    exp.sortSpectra()  
    # !!!

    min_mz: float = sys.maxsize
    max_mz: float = 0
    mslevelcounts: Dict[int,int] = defaultdict(int)

    spectrum_acquisition_metrics_MS1: Dict[str,List[Any]] = defaultdict(list)
    spectrum_acquisition_metrics_MS2: Dict[str,List[Any]] = defaultdict(list)
    spectrum_topn: Dict[str,List[Any]] = defaultdict(list)
    tandem_spectrum_metrics_MS2: Dict[str,List[Any]] = defaultdict(list)
    trap_metrics_MS1: Dict[str,List[Any]] = defaultdict(list)
    trap_metrics_MS2: Dict[str,List[Any]] = defaultdict(list)

    #ActivationMethod look-up dict
    ams = {getattr(ActivationMethod,i): i for i in dir(ActivationMethod) if type(getattr(ActivationMethod,i))==int }

    intens_sum: int = 0
    last_surveyscan_index:int = 0
    for spin, spec in enumerate(exp):
        mslevelcounts[spec.getMSLevel()] += 1
        
        iontraptime = getTrapTime(spec)
        intens_max = spec.get_peaks()[1].max()
        intens_min = spec.get_peaks()[1].min()
        intens_sum = spec.get_peaks()[1].sum()

        if spec.getMSLevel() == 1:
            last_surveyscan_index = spin
            last_surveyscan_intensity = intens_sum
            last_surveyscan_max = intens_max

            spectrum_acquisition_metrics_MS1['RT'].append(spec.getRT())
            spectrum_acquisition_metrics_MS1['SN'].append(getSN_medianmethod(spec))
            spectrum_acquisition_metrics_MS1['peakcount'].append(spec.size())
            spectrum_acquisition_metrics_MS1['int'].append(intens_sum.item())  # .item() for dtype to pytype

            trap_metrics_MS1['RT'].append(spec.getRT())
            trap_metrics_MS1['iontraptime'].append(iontraptime)
            
        if (spec.getMSLevel() == 2):
            if (spec.getPrecursors()[0].getMZ() < min_mz):
                min_mz = spec.getPrecursors()[0].getMZ()
            if (spec.getPrecursors()[0].getMZ() > max_mz):
                max_mz = spec.getPrecursors()[0].getMZ()

            spectrum_acquisition_metrics_MS2['RT'].append(spec.getRT())
            spectrum_acquisition_metrics_MS2['SN'].append(getSN_medianmethod(spec))
            spectrum_acquisition_metrics_MS2['peakcount'].append(spec.size())
            spectrum_acquisition_metrics_MS2['int'].append(intens_sum.item())  # .item() for dtype to pytype
            spectrum_acquisition_metrics_MS2['native_id'].append(spec_native_id(spec))

            trap_metrics_MS2['RT'].append(spec.getRT())
            trap_metrics_MS2['iontraptime'].append(iontraptime)
            trap_metrics_MS2['activation_method'].append(ams.get(list(spec.getPrecursors()[0].getActivationMethods())[0],'unknown'))
            trap_metrics_MS2['activation_energy'].append(spec.getPrecursors()[0].getMetaValue('collision energy'))

            rank = spin - last_surveyscan_index
            spectrum_topn['RT'].append(spec.getRT())
            spectrum_topn['rank'].append(rank)

            # TODO what to do wit dat???
            # precursor_tab['surveyscan_intensity'].append(intens_sum)

            precursor_index = np.searchsorted(exp[last_surveyscan_index].get_peaks()[0], [exp[spin].getPrecursors()[0].getMZ()])[0]
            if precursor_index != np.array(exp[last_surveyscan_index].get_peaks()).shape[1]:
                precursor_err = spec.getPrecursors()[0].getMZ() - np.array(exp[last_surveyscan_index].get_peaks())[:,precursor_index][0]
                precursor_int = np.array(exp[last_surveyscan_index].get_peaks())[:,precursor_index][1]
            else:
                precursor_err = np.nan
                precursor_int = np.nan
                
            tandem_spectrum_metrics_MS2['RT'].append(spec.getRT())
            tandem_spectrum_metrics_MS2['precursor_intensity'].append(precursor_int)  # TODO different from mzid->mzml getPrecursors[0].getIntensity() ? YES, latter one usually zero
            tandem_spectrum_metrics_MS2['precursor_error'].append(precursor_err)
            tandem_spectrum_metrics_MS2['precursor_mz'].append(spec.getPrecursors()[0].getMZ())
            tandem_spectrum_metrics_MS2['precursor_c'].append(spec.getPrecursors()[0].getCharge())

            tandem_spectrum_metrics_MS2['surveyscan_intensity_sum'].append(last_surveyscan_intensity)            
            tandem_spectrum_metrics_MS2['surveyscan_intensity_max'].append(last_surveyscan_max)            

    ## Spectra detail numbers 
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Spectrum acquisition metric values - MS1", value=spectrum_acquisition_metrics_MS1)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Spectrum acquisition metric values - MS2", value=spectrum_acquisition_metrics_MS2)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Spectra topn ranks", value=spectrum_topn)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Tandem spectrum metric values - MS2", value=tandem_spectrum_metrics_MS2)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Trap metric values - MS1", value=trap_metrics_MS1)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Trap metric values - MS2", value=trap_metrics_MS2)
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
                value=[min_mz,max_mz])
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="RT aquisition range", 
                value=[exp[0].getRT(),exp[exp.size()-1].getRT()])
    )

    ## TIC
    tic_tab: Dict[str,List[Any]] = defaultdict(list)
    chroms = exp.getChromatograms()
    for t in chroms:
      if t.getChromatogramType() == oms.ChromatogramSettings.ChromatogramType.TOTAL_ION_CURRENT_CHROMATOGRAM:
        for chro_peak in t:
            tic_tab['RT'].append(chro_peak.getRT())
            tic_tab['int'].append(chro_peak.getIntensity())

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total ion current chromatogram", 
                value=tic_tab)
    )
    # TODO is there a difference between TIC as defined in MS:1000235 and the chromatogram you get from TRP?? In MZML it says its a MS:1000235 (ion current detected in each of a series of mass spectra) but is it?
    # TODO consider collection of spectrum_native_id
    return mzqc.RunQuality(metadata=meta, qualityMetrics=metrics)

def getIDQuality(exp: oms.MSExperiment, pro_ids: List[oms.ProteinIdentification], pep_ids: List[oms.PeptideIdentification], ms2num: int = 0) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    params = pro_ids[0].getSearchParameters()
    # var_mods = params.variable_modifications

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sequence database name", 
                value=pro_ids[0].getSearchParameters().db.decode())
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sequence database version", 
                value=pro_ids[0].getSearchParameters().db_version.decode())
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sequence database taxonomy", 
                value=pro_ids[0].getSearchParameters().taxonomy.decode())
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
    PROTON_MASS_U = 1.00727646677  # Constants::PROTON_MASS_U unavailable

    score_type = pep_ids[0].getScoreType().decode()

    # TODO find a good home for the psi-ms obo in repo
    obo_url = "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"
    with urllib.request.urlopen(obo_url, timeout=10) as obo_in:
        psims = pronto.Ontology(obo_in)

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
        pid = pep_native_id(pepi)
        if pepi.getHits():
            tmp = pepi.getHits()[0]  # TODO apply custom filters and also consider 'pass_threshold'
            identification_scoring_metrics['RT'].append(pepi.getRT())
            identification_scoring_metrics['c'].append(tmp.getCharge())
            identification_scoring_metrics[score_col_name].append(tmp.getScore())

            tw = (tmp.getSequence().getMonoWeight(0,0) + tmp.getCharge() * PROTON_MASS_U) / tmp.getCharge()
            dppm = getMassDifference(tw, pepi.getMZ(), True)      
            identification_accuracy_metrics['RT'].append(pepi.getRT())
            identification_accuracy_metrics['MZ'].append(pepi.getMZ())
            identification_accuracy_metrics['delta_ppm'].append(dppm)
            err = getMassDifference(tw, pepi.getMZ(), False) 
            identification_accuracy_metrics['abs_error'].append(err)

            hydrophobicity_metrics['RT'].append(pepi.getRT())
            hydrophobicity_metrics['gravy'].append(ProtParam.ProteinAnalysis(tmp.getSequence().toUnmodifiedString().decode()).gravy())

            identification_sequence_metrics['RT'].append(pepi.getRT())
            identification_sequence_metrics['peptide'].append(tmp.getSequence().toString().decode().lstrip().rstrip())
            identification_sequence_metrics['target'].append(tmp.getMetaValue('target_decoy').decode().lower() == 'target')
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

def getChargeRatios(identification_scoring_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
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

# TODO target decoy separation and FDR
def getErrorRates(identification_accuracy_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
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

    qs =  np.quantile(identification_accuracy_metrics.value['delta_ppm'], [.25,.5,.75])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15D", 
                value=qs[2]-qs[0])
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm Q1, Q2, Q3", 
                value=list(qs))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm sigma", 
                value=np.std(identification_accuracy_metrics.value['delta_ppm']))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm mean", 
                value=np.mean(identification_accuracy_metrics.value['delta_ppm']))
    )

    np_dppm = np.array(identification_accuracy_metrics.value['delta_ppm'])
    low_out = qs[0]-(1.5*(qs[2]-qs[0]))
    high_out = qs[2]+(1.5*(qs[2]-qs[0]))
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Delta ppm +/-1.5*IQR outlier", 
                value=np.extract((np_dppm<low_out) | (np_dppm>high_out), np_dppm))
    )

    return metrics

def getSamplingRatios(identification_sequence_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
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

def getESIstability(ion_intensity_metric:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    folds = np.true_divide(ion_intensity_metric.value['int'][:-1],ion_intensity_metric.value['int'][1:])
    jumps = len(np.where(folds > 10)[0])
    falls = len(np.where(folds < 1/10)[0])

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

def getMSdensity(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, start_time: datetime.datetime, ms_level: int) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    rts = [start_time + datetime.timedelta(seconds=i) for i in spectrum_acquisition_metrics_MS.value['RT']]
    freq_pm = [len(list(g)) for k, g in itertools.groupby(rts, key=lambda d: d.minute)]

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="fastest frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=max(freq_pm))
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="slowest frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=min(freq_pm))
    )

    f_q1, f_q2, f_q3 =  np.quantile(freq_pm, [.25,.5,.75])
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Q1, Q2, Q3 of frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=[f_q1, f_q2, f_q3])
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sigma of frequency for MS level {ms_level} collection".format(ms_level=ms_level),
                value=np.std(freq_pm))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Mean of frequency for MS level {ms_level} collection".format(ms_level=ms_level),
                value=np.mean(freq_pm))
    )

    np_frq = np.array(freq_pm)
    low_out = f_q1-(1.5*(f_q3-f_q1))
    high_out = f_q3+(1.5*(f_q3-f_q1))
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Frequency for MS level {ms_level} collection +/-1.5*IQR outlier".format(ms_level=ms_level),
                value=np.extract((np_frq<low_out) | (np_frq>high_out), np_frq))
    )

    np_peakdens = np.array(spectrum_acquisition_metrics_MS.value['peakcount'])
    p_q1, p_q2, p_q3 =  np.quantile(np_peakdens, [.25,.5,.75])
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Q1, Q2, Q3 of peak density for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=[p_q1, p_q2, p_q3])
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sigma of peak density for MS level {ms_level} collection".format(ms_level=ms_level),
                value=np.std(np_peakdens))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Mean of peak density for MS level {ms_level} collection".format(ms_level=ms_level),
                value=np.mean(np_peakdens))
    )

    low_out = p_q1-(1.5*(p_q3-p_q1))
    high_out = p_q3+(1.5*(p_q3-p_q1))
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Peak density for MS level {ms_level} collection +/-1.5*IQR outlier".format(ms_level=ms_level),
                value=np.extract((np_peakdens<low_out) | (np_peakdens>high_out), np_peakdens))
    )

    return metrics

def distributionMetricScores(identification_scoring_metrics:mzqc.QualityMetric, score_type:str) -> List[mzqc.QualityMetric]:
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

# numpy array sort:
# sort 2d array xl according to    row n =  xl[:, xl[n].argsort()]
# sort 2d array xl according to column m =  xl[xl[:,m].argsort()]

# TODO split id-free
def getIdentifiedSignalMetrics(tandem_spectrum_metrics_MS2:mzqc.QualityMetric, identification_accuracy_metrics:mzqc.QualityMetric, 
                                spectrum_acquisition_metrics_MS1:mzqc.QualityMetric, spectrum_acquisition_metrics_MS2:mzqc.QualityMetric, 
                                tic_table:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
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


    # ID-FREE
    # DS-3B reimpl.: median( (surv max / prec int) for bottom 50% of all precursors ) 
    np_prec = np_prec[:,np_prec[4].argsort()]
    # Ratio of MS1 maximum to MS1 value at sampling for bottom 50% of analytes by MS1 maximum intensity (1 = sampled at peak maxima)    
    bottom_sampled_prec = np_prec[:,np_prec[4]<np.median(np_prec[4])]
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median ratio of max survey scan intensity over sampled precursor intensity for the bottom (by MS1 max) half of MS2", 
                value=np.median(bottom_sampled_prec[4] / bottom_sampled_prec[2]))
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

def getCoverageRatios(proteinids: oms.ProteinIdentification, 
                peptideids: List[oms.PeptideIdentification], 
                fasta= Dict[str,SeqRecord.SeqRecord]) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    nup = list()
    for p in proteinids.getHits():
        ac = p.getAccession().decode()
        nup.append(oms.ProteinHit(p)) 
        try:
            nup[-1].setSequence(str(fasta[ac].seq))
        except KeyError:
            warnings.warn("Provided fasta file mismatch with provided identifications, aborting.", Warning)
            return [None]

    proteinids.setHits(nup)
    proteinids.computeCoverage(peptideids)

    acc: List[str] = list()
    cov: List[float] = list()
    plen: List[int] = list()
    td: List[str] = list()
    coverage_tab: Dict[str,List[Any]] = defaultdict(list)
    for p in proteinids.getHits():
        coverage_tab['Accession'].append(p.getAccession().decode())
        coverage_tab['Coverage'].append(p.getCoverage())
        coverage_tab['Length'].append(len(p.getSequence().decode()))
        # TODO figure out decoy string by fasta
        coverage_tab['TD'].append('decoy' if 'decoy' in p.getAccession().decode().lower() else 'target')
    
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Protein coverage", 
                value=coverage_tab)
    )

    return metrics

# TODO matching ends?  #internal 
def matchEnzyme(enzre, pepseq):
    matches = np.array([x.start() if x.start()==x.end() else None for x in enzre.finditer(pepseq)])
    is_matched = False
    is_semi = False
    if 1 in matches and len(pepseq)-1 in matches:
        is_matched = True
        internal_matches = len(matches) -2
    elif 1 in matches or len(pepseq)-1 in matches:
        internal_matches = len(matches) - 1
        is_semi = True
    else:
        internal_matches = len(matches)

    return (2 if is_matched else 1 if is_semi else 0, internal_matches)

# TODO which metrics?
def getEnzymeContaminationMetrics(pepl, prol, force_enzymes = False) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    # include all psm actually does not make much sense to assess the enzyme efficiency
    gre = {prol[0].getSearchParameters().digestion_enzyme.getName().decode(): re.compile(prol[0].getSearchParameters().digestion_enzyme.getRegEx().decode())}
    li: List = list()
    oms.EnzymesDB().getAllNames(li)
    ore = {e.decode(): re.compile(oms.EnzymesDB().getEnzyme(e.decode()).getRegEx().decode()) for e in li if e.decode() not in gre and e.decode() != 'no cleavage'}
    
    matches = dict() 
    alt = dict()
    for i,pepi in enumerate(pepl):
        pepi.sort()
        spec_id = pepi.getMetaValue('spectrum_reference').decode() if pepi.metaValueExists('spectrum_reference') else i
        for h in pepi.getHits():
            pepseq = h.getPeptideEvidences()[0].getAABefore() + h.getSequence().toUnmodifiedString().decode() + h.getPeptideEvidences()[0].getAAAfter()
                
            is_matched, internal_matches = matchEnzyme(next(iter(gre.values())), pepseq)
            if is_matched:
                matches[spec_id] = (is_matched, internal_matches)
            if force_enzymes or not is_matched:
                oth_enz_matched = {k: matchEnzyme(v, pepseq) for k,v in ore.items()}
                alt[spec_id] = oth_enz_matched

            # just the first (best) psm
            break
            
    # XMAGHHHEHEQERDHEQEHEHDSLQRP
    # KPNPASMX
    # RPSTNDPTSCCSX

    # enzyme as in file
    # mc with enzyme as in file
    # enzyme as detected
    # mc with enzyme as detected
    # all enzyme matches
    # majority from unique sequences/rank1/allpsm
    # pie plot???

    # contaminations only found with unspecific search!!!

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number of missed cleavages", 
                value=missedcleavages_total)
    )
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Number of missed cleavages in rank1", 
                value=missedcleavages)
    )

    return metrics

# def calcCoverageHelperDatabase():
    # import xml.etree.cElementTree as cet
    # sdbs: List = list()
    # source = "tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.mzid"
    # for event, elem in cet.iterparse(source):
    #     sdb: Dict = dict()
    #     if elem.tag.endswith("SearchDatabase"):
    #         sdb = elem.attrib
    #         sdb['cvs'] = list()
    #         for chil in elem.getchildren():   
    #             if chil.tag.endswith("DatabaseName"):
    #                 for subchil in chil.getchildren():
    #                     # print("#" + subchil)
    #                     if subchil.tag.endswith("cvParam") and 'accession' in subchil.attrib and  subchil.attrib['accession'] == "MS:1001013":
    #                         sdb['databasename'] = subchil.attrib['value']
    #                         # print(subchil.attrib)
    #                         subchil.clear()
    #             elif chil.tag.endswith("cvParam"):
    #                 print(chil.tag)
    #                 sdb['cvs'].append(chil.attrib)
    #                 print(chil.attrib)
    #                 chil.clear()
    #         sdbs.append(sdb)
    #     elif not(elem.tag.endswith("DatabaseName") or elem.tag.endswith("cvParam")):
    #         elem.clear()

def getPeptideLengthMetrics(identification_sequence_metrics:mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    regex_mod = r'(\([^\(]*\))'
    regex_noaa = r'([^A-Za-z])'
    # TODO test this: '.(iTRAQ4plex)M(Oxidation)C(Carbamidomethyl)HNVNR'
    lengths = np.array([len(re.sub(regex_noaa, '', re.sub(regex_mod, '', x))) for x in identification_sequence_metrics.value['peptide'] ])
    
    qs =  np.quantile(lengths, [.25,.5,.75])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identified peptide lengths Q1, Q2, Q3", 
                value=list(qs))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identified peptide lengths sigma", 
                value=np.std(lengths))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identified peptide lengths mean", 
                value=np.mean(lengths))
    )

    low_out = qs[0]-(1.5*(qs[2]-qs[0]))
    high_out = qs[2]+(1.5*(qs[2]-qs[0]))
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identified peptide lengths +/-1.5*IQR outlier", 
                value=np.extract((lengths<low_out) | (lengths>high_out), lengths))
    )
      
    return metrics

def getSNMetrics(spectrum_acquisition_metrics_MS:mzqc.QualityMetric, ms_level: int) -> List[mzqc.QualityMetric]:    
    metrics: List[mzqc.QualityMetric] = list() 
    np_sn = np.array(spectrum_acquisition_metrics_MS.value['SN'])
    
    qs =  np.quantile(np_sn, [.25,.5,.75])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Signal-to-noise ratio Q1, Q2, Q3 for MS level {ms_level} collection".format(ms_level=ms_level),
                value=list(qs))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Signal-to-noise ratio sigma for MS level {ms_level} collection".format(ms_level=ms_level),
                value=np.std(np_sn))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Signal-to-noise ratio mean for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=np.mean(np_sn))
    )

    low_out = qs[0]-(1.5*(qs[2]-qs[0]))
    high_out = qs[2]+(1.5*(qs[2]-qs[0]))
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Signal-to-noise ratio +/-1.5*IQR outlier for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=np.extract((np_sn<low_out) | (np_sn>high_out), np_sn))
    )
      
    return metrics
  
def get_mq_zipped_evidence(url: str) -> Tuple[pandas.DataFrame,pandas.DataFrame]:
    with urllib.request.urlopen(url, timeout=10) as dl:
        with zipfile.ZipFile(io.BytesIO(dl.read())) as z:
            ef = 'evidence.txt'
            pf = 'parameters.txt'
            ld = dict()  # {'ev':'evidence.txt', 'pa':'parameters.txt'}
            dirs = dict()  # {f: {'ev':'evidence.txt', 'pa':'parameters.txt'} }
            for f in z.namelist():
                if(z.getinfo(f).is_dir()):
                    dirs[f] = dict()
                elif f == ef:
                    ld['ev'] = f
                elif f == pf:
                    ld['pa'] = f

            if len(ld) < 2:
                for f in z.namelist():
                    for d in dirs.keys():
                        # exact expected match otherwise oddities like 'SEARCH/._parameters.txt' are picked up
                        if f == d+ef:
                            dirs[d]['ev'] = f
                        elif f == d+pf:
                            dirs[d]['pa'] = f

                dirs = {k:v for k,v in dirs.items() if len(v)>0}
                if len(dirs) > 1:
                    warnings.warn("MQ result zip contains more than one results folder.", Warning)
                elif len(dirs) < 1:
                    warnings.warn("MQ result zip contains no results, even in subfolders.", Warning)
                    return None, None

                ld = next(iter(dirs.values()))
                if len(ld) < 2:
                    warnings.warn("MQ result zip contains no results.", Warning)
                    return None, None

            with z.open(ld['ev']) as e:
                ev = pandas.read_csv(e,sep='\t')
                ev.columns = map(str.lower, ev.columns)
            with z.open(ld['pa']) as p:
                pa = pandas.read_csv(p,sep='\t', dtype={'Parameter':str})
                pa.columns = map(str.lower, pa.columns)
                pa['parameter'] = pa['parameter'].str.lower()
                pa.set_index('parameter', inplace=True)
            return pa,ev

def getMQMetrics(target_raw: str, params: pandas.DataFrame, evidence: pandas.DataFrame, ms2num: int = 0) -> List[mzqc.QualityMetric]: 
    if not target_raw in evidence['raw file'].unique():
        return list()  # TODO warn
    else:
        mq_metrics : List[mzqc.QualityMetric] = list()
        #https://stackoverflow.com/questions/17071871/how-to-select-rows-from-a-dataframe-based-on-column-values
        target_mq = evidence.loc[(evidence['raw file'] == target_raw) & (evidence['ms/ms scan number'].notnull())]

        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Sequence database name", 
                    value=params.loc['fasta file']['value'])
        )

        proteins = len(target_mq['leading proteins'].unique()) 
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Total number of identified proteins", 
                    value=proteins)
        )

        # #     name="Total number of PSM",   # NA
        # metrics.append(
        #     mzqc.QualityMetric(cvRef="QC", 
        #             accession="QC:0000000", 
        #             name="Total number of PSM", 
        #             value=psm_count)
        # )

        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Total number of identified peptide spectra", 
                    value=len(target_mq))
        )

        peptides = len(target_mq['sequence'].unique())
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Total number identified unique peptide sequences", 
                    value=peptides)
        )

        score_type = "Andromeda:score"
        # TODO find a good home for the psi-ms obo in repo
        obo_url = "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"
        with urllib.request.urlopen(obo_url, timeout=10) as obo_in:
            psims = pronto.Ontology(obo_in)

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


        identification_scoring_metrics = target_mq[['retention time','charge','score']].rename(columns={'retention time':'RT','charge': 'c','score':score_type}).to_dict(orient='list')
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Identification scoring metric values", 
                    value=identification_scoring_metrics)
        )

        # TODO comparison column with qccalculator dppm values
        # TODO RT/native id?
        identification_accuracy_metrics = target_mq[['ms/ms m/z','mass error [ppm]','uncalibrated mass error [da]']]\
            .rename(columns={'ms/ms m/z': 'MZ','mass error [ppm]':'delta_ppm','uncalibrated mass error [da]':'abs_error'})
        identification_accuracy_metrics['abs_error'] = identification_accuracy_metrics['abs_error'].abs()
        identification_accuracy_metrics = identification_accuracy_metrics.to_dict(orient='list')    
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Identifications accuracy metric values", 
                    value=identification_accuracy_metrics)
        )

        hydrophobicity_metrics = target_mq[['retention time','sequence']].rename(columns={'retention time':'RT','sequence':'peptide'})
        hydrophobicity_metrics['gravy'] = hydrophobicity_metrics['peptide'].apply(lambda x: ProtParam.ProteinAnalysis(x).gravy())
        hydrophobicity_metrics = hydrophobicity_metrics[['RT','gravy']].to_dict(orient='list')
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Hydrophobicity metric values", 
                    value=hydrophobicity_metrics)
        )

        # TODO target/decoy info available??
        identification_sequence_metrics = target_mq[['sequence','retention time','ms/ms scan number']].rename(columns={'sequence':'peptide','retention time':'RT','ms/ms scan number':'native_id'}).to_dict(orient='list')
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Identifications sequence metric values", 
                    value=identification_sequence_metrics)
        )

        ## simple id metrics 
        mq_metrics.append(
            mzqc.QualityMetric(cvRef="QC", 
                    accession="QC:0000000", 
                    name="Identification to tandem spectra ratio", 
                    value=float(len(target_mq))/float(ms2num))
        )

        return mq_metrics
