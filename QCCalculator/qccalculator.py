import sys
import pyopenms as oms
from MZQC import MZQCFile as mzqc
from typing import List, Dict, Set, Any, Optional, Callable
from collections import defaultdict
from statistics import mean, median, stdev

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

def getBasicQuality(exp: oms.MSExperiment) -> mzqc.RunQuality:
    metrics: List[mzqc.QualityMetric] = list()
    base_name: str = exp.getExperimentalSettings().getSourceFiles()[0].getNameOfFile().decode()
    chksm: str = exp.getExperimentalSettings().getSourceFiles()[0].getChecksum().decode()
    cmpltn: str = exp.getDateTime().get().decode()

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

    exp.sortSpectra()
    min_mz: float = sys.maxsize
    max_mz: float = 0
    mslevelcounts: Dict[int,int] = defaultdict(int)
    
    precursor_tab: Dict[str,List[Any]] = defaultdict(list)

    for ei in exp:
        mslevelcounts[ei.getMSLevel()] += 1
        if (ei.getMSLevel() == 2):
            if (ei.getPrecursors()[0].getMZ() < min_mz):
                min_mz = ei.getPrecursors()[0].getMZ()
            if (ei.getPrecursors()[0].getMZ() > max_mz):
                max_mz = ei.getPrecursors()[0].getMZ()
            precursor_tab['RT'].append(ei.getRT())
            precursor_tab['MZ'].append(ei.getPrecursors()[0].getMZ())
            precursor_tab['c'].append(ei.getPrecursors()[0].getCharge())
            precursor_tab['SN'].append(getSN_medianmethod(ei))
            precursor_tab['peakcount'].append(ei.size())
            # TODO incomplete AcqusitionInfo.pxd
            # it = 'N/A'
            # if not ei.getAcquisitionInfo().empty():
            #     for j in ei.getAcquisitionInfo():
            #         if j.metaValueExists("MS:1000927"):
            #             it = j.getMetaValue("MS:1000927")
            #             break
            # precursor_tab['injectiontime'].append(it)

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

    surveyscan_tab: Dict[str,List[Any]] = defaultdict(list)

    for spec in exp:
        if spec.getMSLevel() == 1:
            intens: int = 0
            for peak in spec:
                intens += peak.getIntensity()
            
            surveyscan_tab['RT'].append(spec.getRT())
            surveyscan_tab['int'].append(intens)
            surveyscan_tab['peakcount'].append(spec.size())
            surveyscan_tab['SN'].append(getSN_medianmethod(spec))

    metrics.append(mzqc.QualityMetric(cvRef="QC", accession="QC:0000000", name="Acquisition table - survey scans", value=surveyscan_tab))

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

    psm_tab: Dict[str,List[Any]] = defaultdict(list)
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
    #varmod???         
    #   for (UInt w = 0; w < var_mods.size(); ++w)
    #   {
    #     at.colTypes.push_back(String(var_mods[w]).substitute(' ', '_'));
    #   }

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
