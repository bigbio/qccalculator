import re
import sys
import datetime 
import itertools
import pyopenms as oms
from MZQC import MZQCFile as mzqc
from typing import List, Dict, Set, Any, Optional, Callable
from collections import defaultdict
from statistics import mean, median, stdev
import numpy as np
from Bio.SeqUtils import ProtParam
from Bio import SeqIO, SeqRecord

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
    
    if (spec.size() % 2 == 0):
        median = (spec[spec.size() // 2].getIntensity() + spec[(spec.size() // 2) +1].getIntensity()) / 2
    else:
        median = spec[spec.size() // 2].getIntensity()
    
    maxi = spec[spec.size()-1].getIntensity()
    if (not norm):
        return maxi / median

    sign_int: float = 0.0
    nois_int: float = 0.0
    sign_cnt: int = 0
    nois_cnt: int  = 0
    for pt in spec:
        if (pt.getIntensity() <= median):
            nois_cnt +=1
            nois_int += pt.getIntensity()
        else:
            sign_cnt +=1 
            sign_int += pt.getIntensity()
      
    return (sign_int / sign_cnt) / (nois_int / nois_cnt)

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

def getBasicQuality(exp: oms.MSExperiment) -> mzqc.RunQuality:
    metrics: List[mzqc.QualityMetric] = list()
    base_name: str = exp.getExperimentalSettings().getSourceFiles()[0].getNameOfFile().decode()
    chksm: str = exp.getExperimentalSettings().getSourceFiles()[0].getChecksum().decode()
    cmpltn: str = exp.getDateTime().get().decode()
    # strt:datetime.datetime = datetime.datetime.strptime(cmpltn, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(seconds=exp.getChromatograms()[0][exp.getChromatograms()[0].size()-1].getRT()*60)

    meta: mzqc.MetaDataParameters = mzqc.MetaDataParameters(
        inputFiles=[
            mzqc.InputFile(name=base_name,location="file:///dev/null", 
                        fileFormat=mzqc.CvParameter("MS", "MS:1000584", "mzML format"), 
                        fileProperties=[
                            mzqc.CvParameter(cvRef="MS", 
                                accession="MS:1000747", 
                                name="completion time", 
                                value=cmpltn
                            ),
                            mzqc.CvParameter(cvRef="MS", 
                                accession="MS:1000569", 
                                name="SHA-1", 
                                value=chksm
                            ),
                            mzqc.CvParameter(cvRef="MS", 
                                accession="MS:1000031", 
                                name="instrument model", 
                                value=exp.getInstrument().getName().decode()
                            )
                        ]
            )
        ], 
        analysisSoftware=[
            mzqc.AnalysisSoftware(cvRef="MS", accession="MS:1000752", name="TOPP software", version="2.4", uri="openms.de")
        ]
    )

    exp.sortSpectra()  # this is mighty important to sort by RT
    min_mz: float = sys.maxsize
    max_mz: float = 0
    mslevelcounts: Dict[int,int] = defaultdict(int)
    
    precursor_tab: Dict[str,List[Any]] = defaultdict(list)  # TODO consider collection of spectrum native_id
    surveyscan_tab: Dict[str,List[Any]] = defaultdict(list)

    intens: int = 0
    # surv: oms.MSSpectrum = None
    last_surveyscan_index = 0
    for spin, spec in enumerate(exp):
        mslevelcounts[spec.getMSLevel()] += 1
        
        iontraptime = getTrapTime(spec)

        if spec.getMSLevel() == 1:
            intens = spec.get_peaks()[1].sum()
            last_surveyscan_index = spin
            intens_max = spec.get_peaks()[1].max()
            surveyscan_tab['RT'].append(spec.getRT())
            surveyscan_tab['int'].append(intens)
            surveyscan_tab['peakcount'].append(spec.size())
            surveyscan_tab['SN'].append(getSN_medianmethod(spec))
            surveyscan_tab['iontraptime'].append(iontraptime)
            
        if (spec.getMSLevel() == 2):
            if (spec.getPrecursors()[0].getMZ() < min_mz):
                min_mz = spec.getPrecursors()[0].getMZ()
            if (spec.getPrecursors()[0].getMZ() > max_mz):
                max_mz = spec.getPrecursors()[0].getMZ()
            precursor_tab['RT'].append(spec.getRT())
            precursor_tab['MZ'].append(spec.getPrecursors()[0].getMZ())
            precursor_tab['c'].append(spec.getPrecursors()[0].getCharge())
            precursor_tab['SN'].append(getSN_medianmethod(spec))
            precursor_tab['peakcount'].append(spec.size())
            precursor_tab['peaksum'].append(spec.get_peaks()[1].sum())
            precursor_tab['surveyscan_intensity'].append(intens)
            precursor_tab['surveyscan_max'].append(intens_max)
            precursor_tab['iontraptime'].append(iontraptime)
            
            prec_ind = np.searchsorted(exp[last_surveyscan_index].get_peaks()[0], [exp[spin].getPrecursors()[0].getMZ()])[0]
            pperr = spec.getPrecursors()[0].getMZ() - np.array(exp[last_surveyscan_index].get_peaks())[:,prec_ind][0]

            precursor_tab['precursor_intensity'].append(np.array(exp[last_surveyscan_index].get_peaks())[:,prec_ind][1])
            precursor_tab['precursor_error'].append(pperr)
            
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Acquisition table - survey scans", value=surveyscan_tab)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Acquisition table - precursors", value=precursor_tab)
    )

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

    tic_tab: Dict[str,List[Any]] = defaultdict(list)
    chroms = exp.getChromatograms()
    for t in chroms:
      if t.getChromatogramType() == oms.ChromatogramSettings.ChromatogramType.TOTAL_ION_CURRENT_CHROMATOGRAM:
        for chro_peak in t:
            tic_tab['RT'].append(chro_peak.getRT() * 60)
            tic_tab['int'].append(chro_peak.getIntensity())

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="TIC", 
                value=tic_tab)
    )

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

    for it in peptides_allhits:
        for c in it:
            if (c == 'K' or c == 'R'):
                missedcleavages_total += 1

    for it in peptides:
        for c in it:
            if (c == 'K' or c == 'R'):
                missedcleavages += 1

    for proid in pro_ids:
        protein_evidence_count += len(proid.getHits())
        for p in proid.getHits():
            proteins.add(p.getAccession())

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
                name="Total number of peptides", 
                value=spectrum_count)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Total number of uniquely identified peptides", 
                value=len(peptides))
    )

    psm_tab: Dict[str,List[Any]] = defaultdict(list)  # TODO consider collection of spectrum native_id
    gravy_tab: Dict[str,List[Any]] = defaultdict(list)
    PROTON_MASS_U = 1.00727646677  # Constants::PROTON_MASS_U  unavailable

    for pepi in pep_ids:
        if pepi.getHits():
            tmp = pepi.getHits()[0]
            psm_tab['RT'].append(pepi.getRT())
            psm_tab['MZ'].append(pepi.getMZ())
            psm_tab['c'].append(tmp.getCharge())
            psm_tab['score'].append(tmp.getScore())
            tw = (tmp.getSequence().getMonoWeight(0,0) + tmp.getCharge() * PROTON_MASS_U) / tmp.getCharge()
            dppm = getMassDifference(tw, pepi.getMZ(), True)
            psm_tab['delta_ppm'].append(dppm)
            psm_tab['peptide_sequence'].append(tmp.getSequence().toString().decode().lstrip().rstrip())
            psm_tab['theoretical_weight'].append(tw)

            gravy_tab['RT'].append(pepi.getRT())
            gravy_tab['RT'].append(ProtParam.ProteinAnalysis(tmp.getSequence().toUnmodifiedString().decode()).gravy())
    #varmod???         
    #   for (UInt w = 0; w < var_mods.size(); ++w)
    #   {
    #     at.colTypes.push_back(String(var_mods[w]).substitute(' ', '_'));
    #   }

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Gravy table", 
                value=gravy_tab)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Identifications table", 
                value=psm_tab)
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Mean delta ppm", 
                value=mean(psm_tab['delta_ppm']))
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median delta ppm", 
                value=median(psm_tab['delta_ppm']))
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="stdev delta ppm", 
                value=stdev(psm_tab['delta_ppm']))
    )

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="id ratio", 
                value=float(len(pep_ids))/float(ms2num))
    )

    return metrics

def getChargeRatios(prec_table) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    # IS3A <- c1n / c2n
    # IS3B <- c3n / c2n
    # IS3C <- c4n / c2n
    unique_charges, charge_freq = np.unique(prec_table.value['c'], return_counts=True)
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
                value=np.median(prec_table.value['c']))
    )

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="IS3Y", 
                value=np.mean(prec_table.value['c']))
    )

    return metrics

def getErrorRates(psm_table) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
 
    err = np.subtract(psm_table.value['MZ'],psm_table.value['theoretical_weight'])

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15A", 
                value= np.median(err) )
    )
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15B", 
                value=np.mean(np.abs(err)) )
    )


    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15C", 
                value=np.median(psm_table.value['delta_ppm']) )
    )

    qs =  np.quantile(psm_table.value['delta_ppm'], [.25,.5,.75])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="MS15D", 
                value=qs[2]-qs[0])
    )

    return metrics

def getSamplingRatios(psm_table) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    sample_rate, sample_rate_counts = np.unique(np.unique(psm_table.value['peptide_sequence'], return_counts=True)[1], return_counts=True)
    explicit_rate_counts = np.zeros( np.max(sample_rate) ) 
    explicit_rate = np.arange(1, np.max(sample_rate)+1)
    explicit_indices = np.where(np.isin(explicit_rate,sample_rate))
    np.put(explicit_rate,explicit_indices,sample_rate)
    np.put(explicit_rate_counts,explicit_indices,sample_rate_counts)
    
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Sampling frequencies", 
                value={'sampling rate': list(explicit_rate), 'frequencies': list(explicit_rate_counts)})
    )
      
    return metrics

def getESIstability(surv_table) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

            # surveyscan_tab['RT'].append(spec.getRT())
            # surveyscan_tab['int'].append(intens)
            # surveyscan_tab['peakcount'].append(spec.size())
            # surveyscan_tab['SN'].append(getSN_medianmethod(spec))

    folds = np.true_divide(surv_table.value['int'][:-1],surv_table.value['int'][1:])
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

def getMSdensity(table: mzqc.QualityMetric, start_time: datetime.datetime, ms_level: int) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    rts = [start_time + datetime.timedelta(seconds=i) for i in table.value['RT']]
    freq_pm = [len(list(g)) for k, g in itertools.groupby(rts, key=lambda d: d.minute)]

    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="fastest frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=max(freq_pm))
    )

    qs =  np.quantile(table.value['peakcount'], [.25,.5,.75])
    metrics.append(
        mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="q1-q3 of frequency for MS level {ms_level} collection".format(ms_level=ms_level), 
                value=max(freq_pm))
    )

    return metrics

def distributionMetricScores(psm_table) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    
    qs =  np.quantile(psm_table.value['score'], [.25,.5,.75])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="scores Q1, Q2, Q3", 
                value=list(qs))
    )

    qs =  np.quantile(psm_table.value['score'], [.25,.5,.75])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="score sigma", 
                value=np.std(psm_table.value['score']))
    )

    return metrics

# numpy array sort:
# sort 2d array xl according to    row n =  xl[:, xl[n].argsort()]
# sort 2d array xl according to column m =  xl[xl[:,m].argsort()]

def getIdentifiedSignalMetrics(prec_table, psm_table, tic_table, surv_table) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 

    # Fraction of total MS2 scans identified in the first quartile of peptides sorted by MS1 maximum intensity
    ps = np.quantile(prec_table.value['precursor_intensity'], [.05,.25,.5,.75,.95])
    np_prec = np.array([prec_table.value['RT'], 
                        prec_table.value['MZ'], 
                        prec_table.value['precursor_intensity'], 
                        prec_table.value['surveyscan_intensity'], 
                        prec_table.value['surveyscan_max']])


    np_prec = np_prec[:,np_prec[2].argsort()]
    # which precursors are in Q1 of precursor_intensity???
    q1int_slice_prec = np_prec[:,np_prec[2]>ps[3]]
    # which of these are identified???
    idcoord = np.array([psm_table.value['RT'],psm_table.value['MZ']])  # TODO make sure intersection is round-proof
    idinter = np.intersect1d(q1int_slice_prec[1],idcoord[1], assume_unique=False, return_indices=True)
    x = q1int_slice_prec[:,idinter[1]]
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Fraction of identified MS2 in the Q1 of precursor intensity.", 
                value=len(x[1])/len(q1int_slice_prec[1]) )
    )


    # Ratio of 95th over 5th percentile MS1 maximum intensity values for identified peptides (approximates dynamic range of signal)
    allidinter = np.intersect1d(np_prec[1],idcoord[1], assume_unique=False, return_indices=True)
    y = np_prec[:,allidinter[1]]
    p05, p95 = np.quantile(y[2], [.05, .95])
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Ratio of 95th over 5th percentile of precursor intensity for identified peptides", 
                value=p95 / p05 )
    )


    # Median maximum MS1 value for identified peptides
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median precursor intensity of identified MS2", 
                value=np.median(y[2]) )
    )


    # MS1-2A Median SN for MS1 spectra in RT range in which half of peptides are identified	
    msn = np.array([surv_table.value['RT'], surv_table.value['SN']])	
    w = msn[:, msn[0]<np.quantile(y[0], [.5])[0] ]
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median SN for MS1 spectra in RT range in which half of peptides are identified", 
                value=np.median(w[1]) )
    )


    # Median TIC value of RT range in which half of peptides are identified	
    tt = np.array([tic_table.value['RT'], tic_table.value['int']])
    z = tt[:, tt[0]<np.quantile(y[0], [.5])[0] ]
    
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median TIC value of RT range in which half of peptides are identified", 
                value=np.median(z[1]) )
    )

    # Ratio of MS1 maximum to MS1 value at sampling for median decile of peptides by MS1 maximum intensity (1 = sampled at peak maxima)
    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median ratio of max survey scan intensity over sampled precursor intensity for peptides identified", 
                value=np.median(y[4] / y[2]))
    )

            
    # Ratio of MS1 maximum to MS1 value at sampling for bottom 50% of peptides by MS1 maximum intensity (1 = sampled at peak maxima)    
    bottom_sampled_prec = np_prec[:,np_prec[4]<np.median(np_prec[4])]

    metrics.append(mzqc.QualityMetric(cvRef="QC", 
                accession="QC:0000000", 
                name="Median ratio of max survey scan intensity over sampled precursor intensity for the bottom (by MS1 max) half of MS2", 
                value=np.median(bottom_sampled_prec[4] / bottom_sampled_prec[2]))
    )

    return metrics

def calcCoverageHelperDatabase():
    import xml.etree.cElementTree as cet
    sdbs: List = list()
    source = "tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.mzid"
    for event, elem in cet.iterparse(source):
        sdb: Dict = dict()
        if elem.tag.endswith("SearchDatabase"):
            sdb = elem.attrib
            sdb['cvs'] = list()
            for chil in elem.getchildren():   
                if chil.tag.endswith("DatabaseName"):
                    for subchil in chil.getchildren():
                        # print("#" + subchil)
                        if subchil.tag.endswith("cvParam") and 'accession' in subchil.attrib and  subchil.attrib['accession'] == "MS:1001013":
                            sdb['databasename'] = subchil.attrib['value']
                            # print(subchil.attrib)
                            subchil.clear()
                elif chil.tag.endswith("cvParam"):
                    print(chil.tag)
                    sdb['cvs'].append(chil.attrib)
                    print(chil.attrib)
                    chil.clear()
            sdbs.append(sdb)
        elif not(elem.tag.endswith("DatabaseName") or elem.tag.endswith("cvParam")):
            elem.clear()

# TODO make it a metric (is plot now)
def getCoverageRatios(proteinids: oms.ProteinIdentification, 
                peptideids: List[oms.PeptideIdentification], 
                fasta= Dict[str,SeqRecord.SeqRecord]) -> List[mzqc.QualityMetric]:
    metrics: List[mzqc.QualityMetric] = list() 
    # coverage: sum(len(pep)) / len(prot)  -  sum(sum(len(pep))) / len(prot)
    # with open("tests/iPRG2015_decoy.fasta","r") as f:
    #     fasta = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    nup = list()
    for p in proteinids.getHits():
        ac = p.getAccession().decode()
        nup.append(oms.ProteinHit(p)) 
        nup[-1].setSequence(str(fasta[ac].seq))

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


# matching ends?  #internal 
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

    return metrics

