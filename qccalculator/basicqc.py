import datetime
import sys
import logging
import configparser
from os.path import basename
from itertools import chain
from collections import defaultdict

from typing import Any, Dict, List, Tuple, Set

import pandas as pd
import numpy as np
import pyopenms as oms
from pyopenms import ActivationMethod
from mzqc import MZQCFile as mzqc
from Bio.SeqUtils import ProtParam

from qccalculator import utils
from qccalculator import noiseqc

"""
Basic set of methods to get quality metric calculations from peak and identification files under way
"""


def getMetricSourceFramesCommon(exp: oms.MSExperiment, config: configparser.ConfigParser, verbose: bool=False) -> mzqc.RunQuality:
    """
    This function calculates the basic QualityMetrics and collects meta information from a mass spectrometry peak file and creates the related RunQuality object.

    Calculated basic QC metrics and proto-metrics necessary to calculate more elaborate QC metrics with additional data (e.g. ID).

    Parameters
    ----------
    exp : oms.MSExperiment
        The mass spectrometry peak file to calculate metrics from
    config : configparser.ConfigParser
        The config parser object from which peak matching tolerance values will be sourced. May differ for spectrum types (see spectrum metadata 'filterstring').
    verbose : bool, optional
        switches on verbose logging, by default False

    Returns
    -------
    mzqc.RunQuality
        A RunQuality object containing the list of metrics calculated and metadata collected, ready for integration into a mzQC file object. 
    """    
    metrics: List[mzqc.QualityMetric] = list()

    # Collect the metadata information to the input file
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
    strt:datetime.datetime = datetime.datetime.strptime(cmpltn, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(seconds=exp.getChromatograms()[0][exp.getChromatograms()[0].size()-1].getRT()*60)

    # Create the metadata section to the input file for the mzQC
    meta: mzqc.MetaDataParameters = mzqc.MetaDataParameters(
        inputFiles=[
            mzqc.InputFile(name=base_name,location=input_loc, 
                        fileFormat=mzqc.CvParameter("MS", "MS:1000584", "mzML format"), 
                        fileProperties=[
                            mzqc.CvParameter( 
                                accession="MS:1000747", 
                                name="completion time", 
                                value=cmpltn
                            ),
                            mzqc.CvParameter( 
                                accession="MS:1000569", 
                                name="SHA-256", 
                                value=chksm
                            ),
                            mzqc.CvParameter( 
                                accession="MS:1000031", 
                                name="instrument model", 
                                value=exp.getInstrument().getName()
                            ),
                            mzqc.CvParameter( 
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
            mzqc.AnalysisSoftware(accession="MS:1000752", name="TOPP software", version=oms.__version__, uri="openms.de")
        ]
    )

    # this is mighty important to sort by RT
    exp.sortSpectra()  

    # RT, SN, peakcount, int_sum (,native_id) (,traptime)
    spectrum_acquisition_metrics_MS1: Dict[str,List[Any]] = defaultdict(list)
    # RT, SN, peakscount, int_sum (,native_id) (,traptime) (,rank) 
    spectrum_acquisition_metrics_MS2: Dict[str,List[Any]] = defaultdict(list)
    # RT, int
    tic_tab: Dict[str,List[Any]] = defaultdict(list)

    # RT, precursor_int, precursor_error, precursor_mz, precursor_c, precursor_int_dev (,native_id)
    tandem_spectrum_metrics: Dict[str,List[Any]] = defaultdict(list)
    # RT, iso_target, iso_lower/upper, frag_energy, frag_method, surveyscan_int_sum, surveyscan_int_max, peakcount_iso, int_sum_iso, int_target_iso, iso_peaks_in_ms2 (,native_id)
    isolation_window_metrics: Dict[str,List[Any]] = defaultdict(list)
    
    # ActivationMethod look-up dict
    ams = {getattr(ActivationMethod,i): i for i in dir(ActivationMethod) if type(getattr(ActivationMethod,i))==int }
    
    min_mz: float = sys.maxsize
    max_mz: float = 0
    mslevelcounts: Dict[int,int] = defaultdict(int)

    intens_sum: np.float = 0
    last_surveyscan_index:int = 0
    last_surveyscan_intensity = intens_sum
    for spin, spec in enumerate(exp):
        mslevelcounts[spec.getMSLevel()] += 1
        
        nid = utils.getSpectrumNativeID(spec)
        iontraptime = utils.getTrapTime(spec)
        intens_max = spec.get_peaks()[1].max()
        intens_min = spec.get_peaks()[1].min()
        intens_sum = spec.get_peaks()[1].sum()

        if spec.getMSLevel() == 1:
            last_surveyscan_index = spin
            last_surveyscan_intensity = intens_sum
            last_surveyscan_max = intens_max

            spectrum_acquisition_metrics_MS1['RT'].append(spec.getRT())
            spectrum_acquisition_metrics_MS1['native_id'].append(nid)
            spectrum_acquisition_metrics_MS1['SN'].append(utils.getSNMedianMethod(spec))
            spectrum_acquisition_metrics_MS1['peakcount'].append(spec.size())
            spectrum_acquisition_metrics_MS1['int_sum'].append(intens_sum.item())  # .item() for dtype to pytype not necessary anymore
            spectrum_acquisition_metrics_MS1['traptime'].append(iontraptime)

            tic_tab['RT'].append(spec.getRT())
            tic_tab['int'].append(intens_sum)
            
        if (spec.getMSLevel() == 2):
            if (spec.getPrecursors()[0].getMZ() < min_mz):
                min_mz = spec.getPrecursors()[0].getMZ()
            if (spec.getPrecursors()[0].getMZ() > max_mz):
                max_mz = spec.getPrecursors()[0].getMZ()

            spectrum_acquisition_metrics_MS2['RT'].append(spec.getRT())
            spectrum_acquisition_metrics_MS2['native_id'].append(nid)
            spectrum_acquisition_metrics_MS2['SN'].append(utils.getSNMedianMethod(spec))
            spectrum_acquisition_metrics_MS2['peakcount'].append(spec.size())
            spectrum_acquisition_metrics_MS2['int_sum'].append(intens_sum.item())  # .item() for dtype to pytype
            spectrum_acquisition_metrics_MS2['traptime'].append(iontraptime)
            rank = spin - last_surveyscan_index
            spectrum_acquisition_metrics_MS2['rank'].append(rank)

            precursor_index = np.searchsorted(exp[last_surveyscan_index].get_peaks()[0], [exp[spin].getPrecursors()[0].getMZ()])[0]
            if precursor_index != np.array(exp[last_surveyscan_index].get_peaks()).shape[1]:
                precursor_err = spec.getPrecursors()[0].getMZ() - np.array(exp[last_surveyscan_index].get_peaks())[:,precursor_index][0]
                precursor_int = np.array(exp[last_surveyscan_index].get_peaks())[:,precursor_index][1]
            else:
                precursor_err = np.nan
                precursor_int = np.nan
                
            tandem_spectrum_metrics['RT'].append(spec.getRT())
            tandem_spectrum_metrics['native_id'].append(nid)
            tandem_spectrum_metrics['precursor_int'].append(precursor_int)  # N.B. also different from mzid getPrecursors[0].getIntensity() ? In mzml this is sually zero and the precursor_int_dev equals precursor_int
            tandem_spectrum_metrics['precursor_error'].append(precursor_err)
            tandem_spectrum_metrics['precursor_mz'].append(spec.getPrecursors()[0].getMZ())
            tandem_spectrum_metrics['precursor_c'].append(spec.getPrecursors()[0].getCharge())  # TODO make sure unknown is represented as 0
            tandem_spectrum_metrics['precursor_int_dev'].append(np.abs(spec.getPrecursors()[0].getIntensity()-precursor_int))

            isolation_window_metrics['RT'].append(spec.getRT())
            isolation_window_metrics['native_id'].append(nid)
            isolation_window_metrics['iso_target'].append(spec.getPrecursors()[0].getMZ())  # https://github.com/OpenMS/OpenMS/blob/d17cc251fd0c4068eb253b03c9fb107897771fdc/src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp#L1992
            isolation_window_metrics['iso_lower'].append(spec.getPrecursors()[0].getIsolationWindowLowerOffset())
            isolation_window_metrics['iso_upper'].append(spec.getPrecursors()[0].getIsolationWindowUpperOffset())
            isolation_window_metrics['surveyscan_int_sum'].append(last_surveyscan_intensity)            
            isolation_window_metrics['surveyscan_int_max'].append(last_surveyscan_max)   
            isolation_window_metrics['frag_method'].append(ams.get(next(iter(spec.getPrecursors()[0].getActivationMethods()), None), np.nan))
            isolation_window_metrics['frag_energy'].append(spec.getPrecursors()[0].getMetaValue('collision energy') if \
                spec.getPrecursors()[0].metaValueExists('collision energy') else np.nan)

            lower = spec.getPrecursors()[0].getMZ() - spec.getPrecursors()[0].getIsolationWindowLowerOffset()
            upper = spec.getPrecursors()[0].getMZ() + spec.getPrecursors()[0].getIsolationWindowUpperOffset()

            s = np.array([(i.getMZ(),i.getIntensity()) for i in exp[last_surveyscan_index]], ndmin = 2)
            s = s[np.where(np.logical_and(s[:, 0]>=lower, s[:, 0]<=upper))[0]]
            isolation_window_metrics['peakcount_iso'].append(np.shape(s)[0]) 
            
            # calc intensity ratio of peaks within isolation window
            int_sort_desc = np.flip(np.argsort(s[:,1]))
            if np.shape(s)[0] > 1:
                isolation_window_metrics['int_ratio_ranked_iso'].append(
                    s[int_sort_desc][:-1,1]/s[int_sort_desc][1:,1][0])  # intensity ratio between top1&2, 2&3, ...
            else:
                isolation_window_metrics['int_ratio_ranked_iso'].append(0)  # bigger is better, though best is 0
            
            isolation_window_metrics['int_sum_iso'].append(np.sum(s[int_sort_desc][:,1]))
            isolation_window_metrics['int_target_iso'].append(spec.getPrecursors()[0].getIntensity())

            if spec.metaValueExists('filter string'):
                ms2_tol = utils.getSpectrumTolerances(spec, config)
            # TODO document used tolerances in a proto/pseudo metric

            # ms2 peaks directly from isolation window?
            unfragmented = np.any([np.isclose(i[0],[x.getMZ() for x in spec], atol=ms2_tol) for i in s])
            isolation_window_metrics['iso_peaks_in_ms2'].append(str(unfragmented))

            # check for all nan columns native id this way, downstream use of tables w/o native_id is obvious to need a differnet matchin method (RT)
            if np.isnan(spectrum_acquisition_metrics_MS1['native_id']).all():  # works only with numbers(+nan) lists!
                spectrum_acquisition_metrics_MS1.pop('native_id')
            if np.isnan(spectrum_acquisition_metrics_MS2['native_id']).all():  # works only with numbers(+nan) lists!
                spectrum_acquisition_metrics_MS2.pop('native_id')
            if np.isnan(tandem_spectrum_metrics['native_id']).all():  # works only with numbers(+nan) lists!
                tandem_spectrum_metrics.pop('native_id')
            if np.isnan(isolation_window_metrics['native_id']).all():  # works only with numbers(+nan) lists!
                isolation_window_metrics.pop('native_id')


    ## internal data tables 
    metrics.append(
        mzqc.QualityMetric(accession="QCc:0000000", 
                name="DateTime of acquisition start", 
                value=strt)
    )
    metrics.append(
        mzqc.QualityMetric(accession="QCc:0000001", 
                name="Spectrum acquisition data - MS1", 
                value=spectrum_acquisition_metrics_MS1)
    )
    metrics.append(
        mzqc.QualityMetric(accession="QCc:0000002", 
                name="Spectrum acquisition data - MS2", 
                value=spectrum_acquisition_metrics_MS2)
    )
    metrics.append(
        mzqc.QualityMetric(accession="QCc:0000003", 
                name="Tandem spectrum data - MS2", 
                value=tandem_spectrum_metrics)
    )
    metrics.append(
        mzqc.QualityMetric(accession="QCc:0000004", 
                name="Isolation window data", 
                value=isolation_window_metrics)
    )

    ## Spectra numbers
    for levels in mslevelcounts.keys():
        if levels == 1: 
            ac = "QC:4000059"
        elif levels == 2: 
            ac = "QC:4000060"
        else: 
            ac = "QC:0000000"
        metrics.append(
            mzqc.QualityMetric(
                    accession=ac, 
                    name="MS{l} count".format(l=str(levels)), 
                    value=mslevelcounts[levels])
        )

    metrics.append(
        mzqc.QualityMetric(accession="QC:4000110", 
                name="Number of chromatograms", 
                value=len(exp.getChromatograms()))
    )

    ## Ranges
    metrics.append(
        mzqc.QualityMetric(accession="QC:4000113", 
                name="MZ aquisition range", 
                value=[min_mz,max_mz])
    )

    metrics.append(
        mzqc.QualityMetric(accession="QC:4000114", 
                name="RT aquisition range", 
                value=[exp[0].getRT(),exp[exp.size()-1].getRT()])
    )

    # TIC
    metrics.append(
        mzqc.QualityMetric(accession="QC:4000069", 
                name="MS1 Total ion current chromatogram", 
                value=tic_tab)
    )

    # Chrom
    chrom_tab: Dict[str,List[Any]] = defaultdict(list)
    chroms = exp.getChromatograms()
    for t in chroms:
      if t.getChromatogramType() == oms.ChromatogramSettings.ChromatogramType.TOTAL_ION_CURRENT_CHROMATOGRAM:
        for chro_peak in t:
            chrom_tab['RT'].append(chro_peak.getRT())
            chrom_tab['int'].append(chro_peak.getIntensity())
        break

    metrics.append(
        mzqc.QualityMetric(accession="QC:0000000", 
                name="Chromatogram", 
                value=chrom_tab)
    )
    # TODO is there a difference between TIC as defined in MS:1000235 and the chromatogram you get from TRP?? In MZML it says its a MS:1000235 (ion current detected in each of a series of mass spectra) but is it?
    # TODO consider collection of spectrum_native_id
    return mzqc.RunQuality(metadata=meta, qualityMetrics=metrics)

# TODO if unfiltered, filter first? If no decoy then what?
def getMetricSourceFramesIdent(pro_ids: List[oms.ProteinIdentification],
                    pep_ids: List[oms.PeptideIdentification], 
                    config: configparser.ConfigParser,
                    id_run_idx: int = 0) -> Tuple[Any,Any,Any]:
    """
    getIDQuality calculates the id-based QualityMetrics from a mass spectrometry peak file and associated 
    identification file.

    Calculated are the id-based QC metrics and proto-metrics necessary to calculate more elaborate QC metrics 
    with even more additional data (e.g. multiple runs).

    Parameters
    ----------
    pro_ids : List[oms.ProteinIdentification]
        List of PyOpenMS ProteinIdentification as from reading a common identification file
    pep_ids : List[oms.PeptideIdentification]
        List of PyOpenMS PeptideIdentification as from reading a common identification file
    ms2num : int, optional
        The total number of tandem spectra as from the id-free metrics, by default 0
    id_run_idx : int, optional
        The id run index used in case multiple id runs are detected, by default 0

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
    metrics: List[mzqc.QualityMetric] = list()
    
    if len(pro_ids) < 1 or id_run_idx > len(pro_ids):
        raise AttributeError("Found NO identification runs at index {}, cannot proceed.".format(str(id_run_idx)))
    elif len(pro_ids) > 1:
        logging.warning("Found multiple identification runs, using index {}.".format(str(id_run_idx))) 

    psims = utils.obtainOntology("psi-ms", config)
    # #pep_ids.sort()

    # TODO there are a number of issues with psims cv subclasses (e.g. MS:1003115 not is_a: MS:1002354 ! PSM-level q-value )
    name_indexed = {psims[x].name: psims[x] for x in psims}
    score_indexed = {x.name: x for x in chain(psims['MS:1001143'].subclasses(), psims['MS:1001153'].subclasses(),
                                              psims['MS:1002347'].subclasses(), psims['MS:1002363'].subclasses())}
    FDRequivalent_indexed =  {x.name: x for x in chain(psims['MS:1002354'].subclasses(),
                                                        [psims['MS:1001868'],psims['MS:1001868'],psims['MS:1001870'],
                                                        psims['MS:1002360'],psims['MS:1002354'],psims['MS:1002355'],
                                                        psims['MS:1003115'],psims['MS:1003116'],psims['MS:1001493'],
                                                        psims['MS:1001491'],psims['MS:1002054'],psims['MS:1002055'],
                                                        psims['MS:1002056'],psims['MS:1001901'],psims['MS:1002265'],
                                                        psims['MS:1001901'],psims['MS:1003098'],psims['MS:1003113']])}

    pep_score_types = set([x.getScoreType() for x in pep_ids])
    pro_score_types = set([x.getScoreType() for x in pro_ids])

    if len(pro_score_types) < 1:
        logging.warning("No protein identification score type set.") 
    if len(pep_score_types) > 1:
        raise AttributeError("Found differing peptide score types assigned to the peptide identifications, cannot proceed.")
        # TODO this can be attempted to be fixed by checking which score types are common and then assign uniformly
    if len(pep_score_types) < 1:
        raise AttributeError("Found NO peptide score types assigned to the peptide identifications, cannot proceed.")
        # TODO this can be attempted to be fixed by checking which score types are common and then assign uniformly
    score_type = next(iter(pep_score_types))

    if score_type in name_indexed:
        if not score_type in score_indexed:
            logging.warning("Score type does not correspond to a score type in the OBO, proceed at own risk.")
            score_col_name = name_indexed[score_type].id
        else:
            score_col_name = score_indexed[score_type].id
    else:
        logging.warn("OBO does not contain any entry matching the identification score, proceed at own risk.")
        score_col_name = score_type

    has_fdr_equiv = False
    if score_type in FDRequivalent_indexed:
        has_fdr_equiv = True

    pep_threshs = set([x.getSignificanceThreshold() for x in pep_ids])
    if len(pep_threshs) > 1:
        logging.warning("Found multiple significance thresholds for peptides, using smallest.")
    pep_sign_thresh = min(pep_threshs)
    pro_threshs = set([x.getSignificanceThreshold() for x in pro_ids])
    if len(pro_threshs) > 1:
        logging.warning("Found multiple significance thresholds for proteins, using smallest.")
    pep_sign_thresh = min(pro_threshs)

    n_id_spec = len(pep_ids)
    n_decoy = 0
    n_psm = 0
    n_no_td_label = 0
    n_r1_decoy = 0
    n_rN = 0
    for pi in pep_ids:
        hits = pi.getHits()
        for i,hi in enumerate(hits):
            if hi.metaValueExists('target_decoy'):
                n_decoy += 1 if 'decoy' in hi.getMetaValue('target_decoy') else 0
                if i == 0 and 'decoy' in hi.getMetaValue('target_decoy'):
                    n_r1_decoy += 1 
            else: 
                n_no_td_label += 1
        n_psm += len(hits) 
        n_rN += 1 if len(hits) > 1 else 0

    if n_no_td_label > 1:
        raise AttributeError("Found PSM without target_decoy label, cannot proceed.")

    filtered = False if n_r1_decoy/n_id_spec > pep_sign_thresh else True
    pruned = False if n_psm > n_id_spec else True

    fdr_filter = max(float(config['fdr.cutoff']['strict']),pep_sign_thresh) 

    if not pruned and not filtered:
        psms = np.ndarray()
        # prune -> peptides
        # filter -> fdr_peptides
    if not filtered and pruned:
        psms = None
        peptides = pd.DataFrame(utils.iterablePSMfromOpenMS(pep_ids),
                        columns=["RT","MZ","c","score","proteinid","sequence","decoy","nativeid"]
                    )
        fdr_peptides = peptides[peptides.score < config['fdr.cutoff']['default']] if has_fdr_equiv else None
    if filtered and pruned:
        psms = None
        peptides = None 
        fdr_peptides = pd.DataFrame(utils.iterablePSMfromOpenMS(pep_ids),
                        columns=["RT","MZ","c","score","proteinid","sequence","decoy","nativeid"]
                    )
    # # RT,MZ,c,score,proteinid,sequence,decoy,nativeid,before,after
    # pd_gen_colname_dtype = pd.DataFrame(utils.iterablePSMfromOpenMS(pep_ids),
    #     columns=["RT","MZ","c","score","proteinid","sequence","decoy","nativeid"],
    #     dtype=np.dtype(object)
    # )



    # params = pro_ids[0].getSearchParameters()
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


    decoy_spectrum_ratio = 0
    # no way to know except for perocolator score
    perco_evidence = None

    spectrum_count: int = 0
    psm_count: int = 0
    runs_count: int = 0
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

    identification_accuracy_metrics: Dict[str, List[Any]] = defaultdict(list)
    identification_scoring_metrics: Dict[str, List[Any]] = defaultdict(list)
    identification_sequence_metrics: Dict[str, List[Any]] = defaultdict(list)
    hydrophobicity_metrics: Dict[str, List[Any]] = defaultdict(list)

    # TODO constants available since 2.5 as oms.Constants.PROTON_MASS_U
    # PROTON_MASS_U = 1.00727646677  # Constants::PROTON_MASS_U unavailable

    for pepi in pep_ids:
        pid = utils.pep_native_id(pepi)
        if pepi.getHits():
            tmp = pepi.getHits()[0]  # TODO apply custom filters and also consider 'pass_threshold'
            identification_scoring_metrics['RT'].append(pepi.getRT())
            identification_scoring_metrics['c'].append(tmp.getCharge())
            identification_scoring_metrics[score_col_name].append(tmp.getScore())

            tw = (tmp.getSequence().getMonoWeight(0, 0) + tmp.getCharge() * oms.Constants.PROTON_MASS_U) / tmp.getCharge()
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
                        value=float(len(pep_ids)) / float(ms2num))
    )

    return psms, peptides, fdr_peptides, metrics
