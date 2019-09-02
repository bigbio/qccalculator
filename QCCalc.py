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


###### rpy2
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.robjects as robjects
from rpy2.robjects import Formula, Environment
from rpy2.robjects.vectors import IntVector, FloatVector
from rpy2.robjects.lib import grid
from rpy2.robjects.packages import importr, data
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects.lib.dplyr import DataFrame
from rpy2.robjects.lib.dplyr import (mutate, group_by, summarize)
from rpy2.rinterface import parse
import warnings
import math, datetime
import tempfile
import base64
import re
import numpy as np
from enum import Enum
import os

class PlotType(Enum):
    PNG = 1
    SVG = 2
    PLOTLY = 3

def handle_plot_format(pp, plot_type: PlotType):
    if plot_type == PlotType.PLOTLY:
        plotly = importr('plotly')
        ppp = plotly.ggplotly(pp)
        htmlwidgets = importr('htmlwidgets')
        with tempfile.NamedTemporaryFile() as t:
                htmlwidgets.saveWidget(ppp, t.name, libdir='lib', selfcontained = False)
                # start stupid fix to get all the recent libs written in the flask lib directory
                htmlwidgets.saveWidget(ppp, 'bof', libdir='lib', selfcontained = False)
                os.remove('bof')
                # end stupid fix
                with open(t.name, "r") as f:
                    s = f.read()
        return s
    else:
        with tempfile.NamedTemporaryFile() as t:
            if plot_type == PlotType.SVG:
                grdevices = importr('grDevices')
                grdevices.svg(file=t.name)
                pp.plot()
                grdevices.dev_off()

                with open(t.name, "r") as f:
                    s = f.read()
            else:
                grdevices = importr('grDevices')
                grdevices.png(file=t.name, width=512, height=512)
                pp.plot()
                grdevices.dev_off()

                with open(t.name, "rb") as fb:
                    s = base64.b64encode(fb.read()).decode()
        return s

def plot_TIC(tic_table, start_time, plot_type=PlotType.PNG):
    d= {'RT': robjects.POSIXct((tuple([start_time + datetime.timedelta(seconds=i) for i in tic_table.value['RT']]))),
        'int': robjects.FloatVector(tuple(tic_table.value['int']))   }
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    c0 = robjects.r('c(0,0)')

    lim_maj=int(max(tic_table.value['RT'])//(60*30))
    lim_min=int(max(tic_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([start_time + datetime.timedelta(seconds=60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([start_time + datetime.timedelta(seconds=60*10* i) for i in range(0,lim_min+1)]))

    axislabels = robjects.StrVector(tuple([(datetime.datetime.fromtimestamp(60*30* i)).strftime("%H:%M") for i in range(0,lim_maj+1)]))

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_line() + \
        ggplot2.aes_string(x='RT', y='int') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = axislabels, expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Intensity", x="Time") + \
        ggplot2.ggtitle("TIC")
    #does not work: date_minor_breaks=scales.date_breaks("5 minutes")
    # scales.date_format("%H:%M")

    # ltb = robjects.r('theme(plot.margin = unit(c(.1,1,.1,.1), "cm"))')
    # pp = pp + ltb

    return handle_plot_format(pp, plot_type)

def plot_SN(table, mslevel=2, svg_plot=False):
    d= {'SN': robjects.FloatVector(tuple(table.value['SN'])) }
    dataf = robjects.DataFrame(d)
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    m = median(table.value['SN'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='SN', y='..density..') + \
        ggplot2.geom_density(alpha=.1, fill="green") + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(SN, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(x="S/N") + \
        ggplot2.ggtitle("S/N distribution in MS{} spectra".format(str(mslevel)))
    
    return handle_plot_format(pp, svg_plot)

def plot_dppm(psm_table, svg_plot=False):
    d= {'deltaPPM': robjects.FloatVector(tuple(psm_table.value['delta_ppm']))   }
    dataf = robjects.DataFrame(d)
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = median(psm_table.value['delta_ppm'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='deltaPPM', y='..density..') + \
        ggplot2.geom_density(alpha=.1, fill="green") + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(deltaPPM, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0, limit=c10) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(x=parse('paste(Delta, "ppm")'), y="Frequency density")  + \
        ggplot2.ggtitle("Mass error distribution")
    
    return handle_plot_format(pp, svg_plot)

def plot_lengths(psm_table, svg_plot=False):
    regex_mod = r'(\([^\(]*\))'
    regex_noaa = r'([^A-Za-z])'
    d= {'PeptideSequence': robjects.StrVector(tuple(psm_table.value['peptide_sequence'])),
        'Length': robjects.IntVector(tuple([len(re.sub(regex_noaa, '', re.sub(regex_mod, '', x))) for x in psm_table.value['peptide_sequence']])) }
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = mean(d['Length'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=1, origin = -0.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='Length') + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='mean(Length, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'mean={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("Length distribution of identified peptide sequences")
        # parse('paste(Delta, "ppm")' does not work in ggplot2.ggtitle

    return handle_plot_format(pp, svg_plot)

def plot_topn(prec_table, surv_table, svg_plot=False):
    h = np.histogram(prec_table.value['RT'], bins=surv_table.value['RT']+[surv_table.value['RT'][-1]+surv_table.value['RT'][-2]])

    d= {'SN': robjects.FloatVector(tuple(surv_table.value['SN'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in surv_table.value['RT']]))), 
        'TopN': robjects.IntVector(tuple(h[0])) }
    dataf = robjects.DataFrame(d)

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = median(d['TopN'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=1, origin = -0.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='TopN') + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(TopN, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(breaks=scales.pretty_breaks(), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("TopN sampling")

    return handle_plot_format(pp, svg_plot)

def plot_topn_sn(prec_table, surv_table, svg_plot=False):
    h = np.histogram(prec_table.value['RT'], bins=surv_table.value['RT']+[surv_table.value['RT'][-1]+surv_table.value['RT'][-2]])
    qs =  np.quantile(surv_table.value['SN'], [.25,.5,.75])
    d= {'SN': robjects.FloatVector(tuple(surv_table.value['SN'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in surv_table.value['RT']]))), 
        'TopN': robjects.IntVector(tuple(h[0])),
        'SN.quartile': robjects.FactorVector((tuple( ['q1' if v < qs[0] else 'q2' if qs[0]<v<qs[1] else 'q3' if qs[1]<v<qs[2] else 'q4' for v in surv_table.value['SN']] )))}
    dataf = robjects.DataFrame(d)

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = median(d['TopN'])
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.aes_string(x='TopN', fill='SN.quartile') + \
        ggplot2.geom_histogram(binwidth=1, origin = -0.5) + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='median(TopN, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'median={}'".format(str(round(m))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(breaks=scales.pretty_breaks(), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("TopN sampling")
        # parse('paste(Delta, "ppm")' does not work in ggplot2.ggtitle

    return handle_plot_format(pp, svg_plot)

def plot_topn_rt(prec_table, surv_table, svg_plot=False):
    h = np.histogram(prec_table.value['RT'], bins=surv_table.value['RT']+[surv_table.value['RT'][-1]+surv_table.value['RT'][-2]])

    d= {'SN': robjects.FloatVector(tuple(surv_table.value['SN'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in surv_table.value['RT']]))), 
        'TopN': robjects.IntVector(tuple(h[0])) }
    dataf = robjects.DataFrame(d)

    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(surv_table.value['RT'])//(60*30))
    lim_min=int(max(surv_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_col(width = 1) + \
        ggplot2.aes_string(x='RT',y='TopN') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = scales.date_format("%H:%M"), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count", x="TopN sampling range")  + \
        ggplot2.ggtitle("TopN utilisation")


    # TODO also plot real histogram with color of SN or target/decoy?
    return handle_plot_format(pp, svg_plot)

def plot_idmap(prec_table, psm_table, svg_plot=False):
    d_psm= {'MZ': robjects.FloatVector(tuple(psm_table.value['MZ'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in psm_table.value['RT']]))), 
        'col': robjects.FactorVector(tuple(["identified"]*len(psm_table.value['MZ']))) }
    dataf_psm = robjects.DataFrame(d_psm)

    d_prc= {'MZ': robjects.FloatVector(tuple(prec_table.value['MZ'])), 
        'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in prec_table.value['RT']]))), 
        'col': robjects.FactorVector(tuple(["recorded"]*len(prec_table.value['MZ']))) }
    dataf_prc = robjects.DataFrame(d_prc)

    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(prec_table.value['RT'])//(60*30))
    lim_min=int(max(prec_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    pp = ggplot2.ggplot(dataf_prc) + \
        ggplot2.geom_point() + \
        ggplot2.aes_string(x='RT',y='MZ', colour='col') + \
        ggplot2.scale_x_datetime(breaks=b_maj, minor_breaks=b_min, labels = scales.date_format("%H:%M"), expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="m/z", x="Time")  + \
        ggplot2.geom_point(data=dataf_psm) + \
        ggplot2.ggtitle("ID map")

    ltb = robjects.r('theme(legend.title=element_blank())')
    pp = pp + ltb

    # TODO also plot real histogram with color of SN or target/decoy?
    return handle_plot_format(pp, svg_plot)

def plot_charge(prec_table, psm_table=None, svg_plot=False):
    d_prc= {'c': robjects.IntVector(tuple(prec_table.value['c'])), 
        'col': robjects.FactorVector(tuple(["recorded"]*len(prec_table.value['c']))) }
    dataf = robjects.DataFrame(d_prc)

    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_bar() + \
        ggplot2.aes_string(x='c',fill='col') + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.scale_x_continuous(breaks=scales.pretty_breaks(), expand=c0) + \
        ggplot2.labs(x="Charge state", y="Count")  + \
        ggplot2.ggtitle("Charge states") + \
        robjects.r('theme(legend.title=element_blank())')

    # TODO N/A handling???

    if psm_table:
        d_prc= {'c': robjects.IntVector(tuple(psm_table.value['c'])), 
            'col': robjects.FactorVector(tuple(["identified"]*len(psm_table.value['c']))) }
        dataf_id = robjects.DataFrame(d_prc)
        pp = pp + ggplot2.geom_bar(data=dataf_id)

    # TODO also plot real histogram with color of SN or target/decoy?
    return handle_plot_format(pp, svg_plot)

def plot_peaknum(table, mslevel=2, svg_plot=False):
    d= {'peakcount': robjects.IntVector(tuple(table.value['peakcount']))}
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')
    m = mean(d['peakcount'])
    binw = round(m/15)+1  #+1 avoids 0 binwidth
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(binwidth=binw, origin = -0.5, colour="black", fill="white") + \
        ggplot2.aes_string(x='peakcount') + \
        ggplot2.geom_vline(ggplot2.aes(xintercept='mean(peakcount, na.rm=TRUE)'), color="red", linetype="dashed", size=1) + \
        ggplot2.geom_text(ggplot2.aes_string(x=str(m), y=rinf, label="'mean={}'".format(str(round(m,2))) ), hjust="left", vjust="top") + \
        ggplot2.scale_x_continuous(expand=c0) + \
        ggplot2.scale_y_continuous(expand=c0) + \
        ggplot2.labs(y="Count")  + \
        ggplot2.ggtitle("Peak count distribution for MS{} spectra".format(str(mslevel)))

    return handle_plot_format(pp, svg_plot)

def plot_intensities(table, svg_plot=False):
    grdevices = importr('grDevices')
    d= {'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in table.value['RT']]))),
        'int': robjects.FloatVector(tuple(table.value['int']))}
    dataf = robjects.DataFrame(d)
    scales = importr('scales')
    cow = importr('cowplot')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    combiaxis1 = robjects.r('theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin=unit(c(1,1,-0.1,1), "cm"))')
    combiaxis2 = robjects.r('theme(axis.text.y = element_blank(), axis.ticks = element_blank(), plot.margin=unit(c(-0.1,1,1,1), "cm"))')

    ht = ggplot2.ggplot(dataf) + \
        ggplot2.geom_histogram(bins=100) + \
        ggplot2.aes_string(x='int') + \
        combiaxis1 + \
        ggplot2.labs(y="Count", x= robjects.r('element_blank()') )  + \
        ggplot2.ggtitle("Intensity distribution")

    bx = ggplot2.ggplot(dataf) + \
        ggplot2.geom_boxplot() + \
        ggplot2.aes_string(x=1, y='int') + \
        ggplot2.coord_flip() + \
        combiaxis2 + \
        ggplot2.labs(x="", y="Intensities")
        
    # pracma = importr('pracma')
    # AUC <- trapz(QCTIC$MS.1000894_.sec.,QCTIC$MS.1000285)
    # auc = np.trapz(table.value['int'], x=table.value['RT'])

    # Qs <- quantile(QCTIC$MS.1000285,prob = c(0.25, 0.5, 0.75))
    # Qs <- data.frame(Qs)
    # qs =  np.quantile(table.value['int'], [.25,.5,.75])

    pp = cow.plot_grid(ht,bx, ncol = 1, align = 'v', axis = 'l', rel_heights = robjects.r('c(1,.25)'))


    # grdevices.png(file="tests/grid_png_func.png", width=512, height=512)
    # c = cow.plot_grid(ht,bx, ncol = 1, align = 'v', axis = 'l', rel_heights = robjects.r('c(1,.25)'))
    # c.plot()
    # grdevices.dev_off()

    return handle_plot_format(pp, svg_plot)

def plot_events(tic_table, surv_table, prec_table, psm_table=None, svg_plot=False):
    datasources = [("Chromatogram", tic_table),("MS1",surv_table),("MS2", prec_table)]
    if psm_table:
        datasources.append(("Identifications",psm_table))
    
    dataf = robjects.DataFrame({})
    for annot, tab in datasources:
        if annot == "Chromatogram":
            #quartiles of total intensity recorded
            tic = np.sum(tab.value['int'])
            vec = (np.cumsum(tab.value['int'])/tic)
            qs =  [.25,.5,.75]
        else:
            vec = [i[0]+1 for i in enumerate(tab.value['RT'])]
            qs =  np.quantile(vec, [.25,.5,.75])

        d_c= {'Time': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in tab.value['RT']]))),
            'Quartile': robjects.FactorVector((tuple( ['q1' if v < qs[0] else 'q2' if qs[0]<v<qs[1] else 'q3' if qs[1]<v<qs[2] else 'q4' for v in vec] )))
        }
        td = (DataFrame(d_c).
                group_by('Quartile').
                summarize(n='n()', mx='max(Time)', mi='min(Time)').
                mutate(group="'"+annot+"'", dt='hms::as.hms(mx-mi)') )
        dataf = dataf.rbind(td)

    dataf = ( DataFrame(dataf).arrange('Quartile') )
    scales = importr('scales')
    chron = importr('chron')
    rinf = robjects.r('Inf')
    c0 = robjects.r('c(0,0)')
    c10 = robjects.r('c(-10,10)')

    # TODO time axis fix and Hz annotations
    b=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10)]))
    lim_maj=int(max(tic_table.value['RT'])//(60*30))
    lim_min=int(max(tic_table.value['RT'])//(60*10))+1
    b_maj=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*30* i) for i in range(0,lim_maj+1)]))
    b_min=robjects.POSIXct(tuple([datetime.datetime.fromtimestamp(60*10* i) for i in range(0,lim_min+1)]))

    # TODO shitty time scale still does not work with stacked bars
    pp = ggplot2.ggplot(dataf) + \
        ggplot2.geom_bar(stat='identity', position = ggplot2.position_fill(reverse = True)) + \
        ggplot2.aes_string(x='group',y='mx',fill='Quartile') + \
        ggplot2.coord_flip() + \
        chron.scale_y_chron(format ="%H:%M:%S") + \
        ggplot2.labs(x="Events", y="Time")  + \
        ggplot2.ggtitle("Quartiles of Chromatographic, MS1, MS2, and identification events over RT") 

    dataf.to_csvfile('tests/events.csv')
    return handle_plot_format(pp, svg_plot)
        # ggplot2.scale_y_datetime(breaks=b_maj, minor_breaks=b_min,labels = scales.date_format("%H:%M"), expand=c0) + \


##### test code
exp = oms.MSExperiment()
oms.MzMLFile().load("tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.trfr.mzML", exp)
rq = getBasicQuality(exp)

cmpltn: str = exp.getDateTime().get().decode()
strt:datetime.datetime = datetime.datetime.strptime(cmpltn, '%Y-%m-%d %H:%M:%S') - datetime.timedelta(seconds=exp.getChromatograms()[0][exp.getChromatograms()[0].size()-1].getRT()*60)

srq_d = mzqc.JsonSerialisable().ToJson(rq)
srq_m = mzqc.JsonSerialisable().ToJson(rq, readability=1)
with open("tests/basicrq_mzenc.json", "w") as f:
    f.write(srq_m)

pepl: List[oms.PeptideIdentification] = list()
prol: List[oms.ProteinIdentification] = list()
oms.MzIdentMLFile().load("tests/cptac_itraq_example.mzid", prol, pepl)
ms2s: Any = -1
for qm in rq.qualityMetrics:
    if qm.name == "Number of MS2 spectra":
        ms2s = qm.value
        break

mid = getIDQuality(exp, prol, pepl, ms2s)
rq.qualityMetrics.extend(mid)

cv1 = mzqc.ControlledVocabulary(ref="QC", name="QC", uri="www.eff.off")
cv2 = mzqc.ControlledVocabulary(ref="MS", name="MS", uri="www.eff.off")

oqc = mzqc.MzQcFile(version="0_0_11", runQualities=[rq], controlledVocabularies=[cv1,cv2])
mzq_m = mzqc.JsonSerialisable().ToJson(oqc, readability=1)
with open("tests/mzenc.json", "w") as f:
    f.write("{ \"mzQC\": " + mzq_m + " }")


svg = plot_TIC(rq.qualityMetrics[6],strt,PlotType.SVG)
with open("tests/tic_svg_func.svg","w") as f:
    f.write(svg)

tic = plot_TIC(rq.qualityMetrics[6],strt)
with open("tests/tic_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(tic.encode()))

sn = plot_SN(rq.qualityMetrics[0])
with open("tests/sn_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(sn.encode()))

dppm = plot_dppm(rq.qualityMetrics[18])
with open("tests/dppm_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(dppm.encode()))

topn = plot_topn(rq.qualityMetrics[0], rq.qualityMetrics[7])
with open("tests/topn_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(topn.encode()))

topnsn = plot_topn_sn(rq.qualityMetrics[0], rq.qualityMetrics[7])
with open("tests/topn_sn_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(topnsn.encode()))

idmap = plot_idmap(rq.qualityMetrics[0], rq.qualityMetrics[18])
with open("tests/idmap_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(idmap.encode()))

charges = plot_charge(rq.qualityMetrics[0], rq.qualityMetrics[18])
with open("tests/charges_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(charges.encode()))

pk = plot_peaknum(rq.qualityMetrics[7], 1)
with open("tests/pk1_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(pk.encode()))

pk = plot_peaknum(rq.qualityMetrics[0], 2)
with open("tests/pk2_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(pk.encode()))

gr = plot_intensities(rq.qualityMetrics[6])
with open("tests/int_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(gr.encode()))

ev = plot_events(rq.qualityMetrics[6], rq.qualityMetrics[7], rq.qualityMetrics[0], rq.qualityMetrics[18])
with open("tests/ev_png_func.png","wb") as fb:
    fb.write(base64.decodebytes(ev.encode()))


tab = rq.qualityMetrics[6]
d= {'RT': robjects.POSIXct((tuple([datetime.datetime.fromtimestamp(i) for i in tab.value['RT']]))),
    'int': robjects.FloatVector(tuple(tab.value['int']))   }
dataf = robjects.DataFrame(d)
ap = ggplot2.ggplot(dataf) + \
        ggplot2.geom_line() + \
        ggplot2.aes_string(x='RT', y='int', color='int')

tad = importr('AnomalyDetection')
anomalies = tad.detect_anoms(dataf, num_obs_per_period=1000)
with open("/tmp/ap.png","wb") as fb:
     fb.write(base64.decodebytes(handle_plot_format(anomalies.plot, False).encode()))

ap=tad.time_decompose(df, frequency = "auto", trend = ".5 hrs").anomalize('remainder').plot_anomaly_decomposition()


bsa = importr('anomalize')
