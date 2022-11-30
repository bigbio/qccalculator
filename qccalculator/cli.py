#!/usr/bin/env python
import os
from os.path import basename
from os.path import join
from os.path import dirname
from os.path import expanduser
from os.path import exists
import click
import logging
from datetime import datetime
import gzip
from typing import List
import configparser
import tempfile
import urllib.request
import pandas

import pyopenms as oms
from click import command
from mzqc import MZQCFile as qc

from qccalculator import basicqc, idfree, idqc, utils
from qccalculator.maxquant_analysis_result import MQanalysisResult

rqs: List[qc.RunQuality] = list()
sqs: List[qc.SetQuality] = list()
out = str()
zp = False

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def print_help():
    """
    Print the help of the tool
    :return:
    """
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()


def global_variables(output, zip):
    """
    Calculate quality metrics for given files.
    Multiple files input is possible (each after a "full/basic" COMMAND).
    All metrics of one QCCalculator execution will be stored in on output file.
    If you need separate mzQC files, please execute separately.
    For more information on the different COMMAND types, try QCCalculator COMMAND --help

    Parameters
    ----------
    output: path to the output file qcml
    zip: add zip to the output file

    Returns
    -------

    """
    logging.info("Received output destination {}".format(output))
    global out
    out = output
    if zip:
        global zp
        zp = True


def mzqc_assembly(rqs, sqs):
    """
    Setup the latest version of qcML ontology file

    Parameters
    ----------
    rqs: Run qualities
    sqs: Set Qualities

    Returns
    -------

    """
    #TODO clean 'QCc' proto metrics
    cv_qc = qc.ControlledVocabulary(name="Proteomics Standards Initiative Quality Control Ontology",
                                    version="0.1.0",
                                    uri="https://github.com/HUPO-PSI/qcML-development/blob/master/cv/v0_1_0/qc-cv.obo")
    cv_ms = qc.ControlledVocabulary(name="Proteomics Standards Initiative Mass Spectrometry Ontology",
                                    version="4.1.7", 
                                    uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo")

    return qc.MzQcFile(version="1.0.0",
                        creationDate=datetime.now().isoformat(),
                        runQualities=rqs,
                        setQualities=sqs,
                        controlledVocabularies=[cv_qc, cv_ms])


# @click.pass_context
def finale():
    """
    Gunzip and compress the output files.
    """
    logging.info("Calculated metrics from {} different input peak files".format(len(rqs)))
    logging.info("Attempting to write results to {}{}".format(out, ".gz" if zp else ""))
    if any(rqs) or any(sqs) or out:
        if zp:
            with gzip.GzipFile(out + '.gz', 'w') as fh:
                fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs), readability=1).encode('utf-8'))
        else:
            with open(out, 'w') as fh:
                fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs), readability=0))
    logging.info("Done. Thank you for choosing QCCalculator!")


def common_options(function):
    function = click.option('--mzml', required=True, type=click.Path(), 
                            help="Path to the mzml to compute the quality metrics")(function)
    function = click.option('--output', type=click.Path(), default="out.mzQC",
                            help="The path and name of the desired output file.")(function)
    function = click.option('--zip', default=False,
                            help="Apply gzip to the output. Appends '.gz' to the target filename and applies pretty formatting inside the mzQC json.",
                            is_flag=True)(function)
    function = click.option('--loglevel', default="WARN", type=click.Choice(["INFO","WARN","DEBUG","MUTE"]),
                            help="Level of detail for log feed. (Valid: INFO/WARN/DEBUG/MUTE)")(function)
    function = click.option('--logfile', type=click.Path(), 
                            help="Log output will be diverted to this file. Log level DEBUG will automatically write to prompt, too.")(function)
    return function


# todo add log options to config
def getconfig():
    home = expanduser("~")
    config_file = join(home, '.qccalculator', 'config.ini')
    if not exists(config_file):
        from qccalculator.configfile import config_content
        qccalc_dir = join(home, '.qccalculator')
        if not exists(qccalc_dir):
            os.makedirs(qccalc_dir)
        with open(config_file, 'w') as f:
            f.write(config_content)
        logging.info("Created configfile: {config_file}".format(config_file=config_file))
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def setlog(loglevel, logfile):
    lev = {"INFO": logging.INFO, "WARN": logging.WARN, "DEBUG": logging.DEBUG}.get(loglevel,logging.WARN)
    han = list()
    if loglevel != "MUTE":
        han.append(logging.StreamHandler())
    if logfile:
        han.append(logging.FileHandler(logfile))
    logging.basicConfig(
        level=lev,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=han
    )
    logging.debug("set logger")


@click.command('basic', short_help='Compute QC metrics for an mzML file', context_settings={'show_default': True})
@common_options
def basic(mzml, output, zip, loglevel, logfile):
    """
    Calculate the basic metrics available from virtually every mzML file.
    """
    if mzml is None or output is None:
        print_help()

    config = getconfig()
    setlog(loglevel, logfile)

    exp = oms.MSExperiment()
    oms.MzMLFile().load(click.format_filename(mzml), exp)
    qms = basicqc.getMetricSourceFramesCommon(exp, config)
    
    # value variables in this code block are usually QC metrics (i.e. value in resp. attribute)
    starttime = next(iter(list(filter(lambda x: (x.accession == "QCc:0000000"), qms))), None)
    spectrum_acquisition_metrics_MS1 = next(iter(list(filter(lambda x: (x.accession == "QCc:0000001"), qms))), None)
    spectrum_acquisition_metrics_MS2 = next(iter(list(filter(lambda x: (x.accession == "QCc:0000002"), qms))), None)
    tandem_spectrum_metrics = next(iter(list(filter(lambda x: (x.accession == "QCc:0000003"), qms))), None)
    isolation_window_metrics = next(iter(list(filter(lambda x: (x.accession == "QCc:0000004"), qms))), None) 
    tic_tab = next(iter(list(filter(lambda x: (x.accession == "QC:4000069"), qms))), None) 
    
    metrics2add = list()

    # todo catch empties!

    if spectrum_acquisition_metrics_MS1:
        if starttime:
            metrics2add.extend(idfree.calcMSFrequency(spectrum_acquisition_metrics_MS1, start_time=starttime.value, ms_level=1))
            metrics2add.extend(idfree.calcMSDensity(spectrum_acquisition_metrics_MS1, start_time=starttime.value, ms_level=1))
        metrics2add.extend(idfree.calcSNMetrics(spectrum_acquisition_metrics_MS1, ms_level=1))
        metrics2add.extend(idfree.calcTrapTimeMetrics(spectrum_acquisition_metrics_MS1, ms_level=1))
        metrics2add.extend(idfree.calcMSSpread(spectrum_acquisition_metrics_MS1, ms_level=1))
    if spectrum_acquisition_metrics_MS2:
        if starttime:
            metrics2add.extend(idfree.calcMSFrequency(spectrum_acquisition_metrics_MS2, start_time=starttime.value, ms_level=2))
            metrics2add.extend(idfree.calcMSDensity(spectrum_acquisition_metrics_MS2, start_time=starttime.value, ms_level=2))
        metrics2add.extend(idfree.calcSNMetrics(spectrum_acquisition_metrics_MS2, ms_level=2))
        metrics2add.extend(idfree.calcTrapTimeMetrics(spectrum_acquisition_metrics_MS2, ms_level=2))
        metrics2add.extend(idfree.calcMSSpread(spectrum_acquisition_metrics_MS2, ms_level=2))
    if tandem_spectrum_metrics:
        metrics2add.extend(idfree.calcAnalysedSignalMetrics(tandem_spectrum_metrics, isolation_window_metrics))
        metrics2add.extend(idfree.calcPrecursorIntensityMetrics(tandem_spectrum_metrics))
        metrics2add.extend(idfree.calcPrecursorChargeMetrics(tandem_spectrum_metrics))
    if tic_tab:
        metrics2add.extend(idfree.calcESIStability(tic_tab))
        metrics2add.extend(idfree.calcTICSpread(tic_tab))
    
    # todo either register QC calc as software cv term or add custom CV to the system
    met = qc.MetaDataParameters(
        analysisSoftware=[qc.AnalysisSoftware(accession="QCc:12345678", name="qccalculator")],
        inputFiles=[qc.InputFile(fileFormat=qc.CvParameter(accession="MS:1000584", name="mzML format"), 
                                location=os.path.dirname(mzml), name=os.path.basename(mzml))],)
    rq = qc.RunQuality(metadata=met , qualityMetrics=metrics2add)
    rqs.append(rq)

    logging.warn("Adding {} calculated metrics.".format(len(metrics2add)))

    global_variables(output, zip)
    finale()


@click.command("full", short_help='Compute the QC metrics for spectra file (mzML) and identification data')
@click.option('--mzid', type=click.Path(exists=True),
                help="If you have a corresponding mzid file you need to pass it, too. Mutually exclusive to idxml.")
@click.option('--idxml', type=click.Path(exists=True),
                help="If you have a corresponding idxml file you need to pass it, too. Mutually exclusive to mzid.")
@common_options
def full(mzid=None, idxml=None, mzml=None, output=None, zip=None, loglevel=None, logfile=None):
    """
    Calculate all possible metrics for these files. These data sources will be included in set metrics.
    """

    if (mzml is None) or (mzid is None and idxml is None):
        print_help()

    config = getconfig()
    setlog(loglevel, logfile)

    exp = oms.MSExperiment()
    oms.MzMLFile().load(click.format_filename(mzml), exp)
    rq = basicqc.getMetricSourceFramesCommon(exp, config)

    if idxml and mzid:
        with click.Context(command) as ctx:
            logging.warn("Sorry, you can only give one id file. Please choose one.")
            click.echo(command.get_help(ctx))
            return
    elif not idxml and not mzid:
        with click.Context(command) as ctx:
            logging.warn("Sorry, you must give one id file in this mode.")
            click.echo(command.get_help(ctx))
            return

    ms2num = 0
    for x in rq.qualityMetrics:
        if x.name == "MS2 count":
            ms2num = x.value

    if ms2num < 1:
        logging.warn("We seem to have found no MS2 spectra which is unlikely to be true since you have also given some identifications. \
                                We continue with symbolic value of 1 for the number of MS2 spectra, \
                                however this means some metrics will invariably be incorrect!\
                                Please make sure, we have the right inputs.")
        ms2num = 1

    pros = list()
    peps = list()
    if mzid:
        oms_id = oms.MzIdentMLFile()
        idf = mzid
    if idxml:
        oms_id = oms.IdXMLFile()
        idf = idxml
    if idf:
        oms_id.load(click.format_filename(idf), pros, peps)
        rq.qualityMetrics.extend(idqc.getIDQuality(pro_ids=pros, pep_ids=peps, ms2num=ms2num, config=config))
    rqs.append(rq)

    global_variables(output, zip)
    finale()


def basics_from_mzml(mzmlpath: str, config: configparser.ConfigParser) -> List[qc.QualityMetric]:
    """
    idfree_from_mzml generates id free source frames (proto-metrics) from one mzML file

    The mzML can be provided via an URI (either http://, https://, or file://).

    Parameters
    ----------
    mzmlpath : str
        The mzML file location.

    Returns
    -------
    List[mzqc.QualityMetric]
        List of (proto-) metrics for the given run.
    """
    exp = oms.MSExperiment()
    sha = None
    if mzmlpath.startswith('http'):
        with tempfile.NamedTemporaryFile() as tf:
            urllib.request.urlretrieve(mzmlpath, tf.name)
            oms.MzMLFile().load(tf.name, exp)
            exp.getExperimentalSettings().setLoadedFilePath(mzmlpath)
            sha = utils.sha256fromfile(tf.name)
    else:
        usepath = utils.remove_protocol_prefix(mzmlpath)
        oms.MzMLFile().load(usepath, exp)
    rq = basicqc.getMetricSourceFramesCommon(exp, config)
    return rq


@click.command("maxquant", short_help="Compute quality metrics from MaxQuant results")
@click.option('--mqpar', type=str, required=True,
                help="The URI to a max quant parameter XML file.")
@click.option('--zipuri', type=str, required=True,
                help="The URI to a max quant output zip file (must contain evidence.txt and parameters.txt).")
@click.option('--eduri', type=str, required=True,
                help="The experimental design in SDRF format.")
@click.option('--mzmlfolder', type=str, default="",
                help="The absolute path or base URL to the mzML files matching the raw files used during MQ analysis.")
@click.option('--output', type=click.Path(), default="out.mzQC",
                help="The path and name of the desired output file.")
@click.option('--zip', default=False, is_flag=True,
                help="Apply gzip to the output. Appends '.gz' to the target filename and applies pretty formatting inside the mzQC json.")
def maxquant(mqpar, zipuri, eduri, mzmlfolder, output, zip, loglevel, logfile):
    """
    Calculate all possible metrics for MaxQuant output. These data sources will be included in set metrics.
    """
    if (zipuri is None or 
        mqpar is None or 
        eduri is None):
        print_help()

    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), 'config.ini'))  # TODO assure access from packaged project
    setlog(loglevel, logfile)

    MQR = MQanalysisResult(mqpar=mqpar, mqzip=zipuri, mqtxt=eduri, peakbase=mzmlfolder)
    for row in MQR.experimentalGrouping.itertuples(index=False):
        rn = row[MQR.experimentalGrouping.columns.get_loc('run name')] 
        mz = row[MQR.experimentalGrouping.columns.get_loc('mzML path')]
        basic = basics_from_mzml(mz, config)
        MQR.run_qualities[rn].qualityMetrics.extend(basic)
    
    for runquality in MQR.run_qualities.values():
        starttime = next(iter(list(filter(lambda x: (x.accession == "QCc:0000000"), runquality.qualityMetrics))), None)
        spectrum_acquisition_metrics_MS1 = next(iter(list(filter(lambda x: (x.accession == "QCc:0000001"), runquality.qualityMetrics))), None)
        spectrum_acquisition_metrics_MS2 = next(iter(list(filter(lambda x: (x.accession == "QCc:0000002"), runquality.qualityMetrics))), None)
        tandem_spectrum_metrics = next(iter(list(filter(lambda x: (x.accession == "QCc:0000003"), runquality.qualityMetrics))), None)
        tic_tab = next(iter(list(filter(lambda x: (x.accession == "QC:4000069"), runquality.qualityMetrics))), None) 
        isolation_window_metrics = next(iter(list(filter(lambda x: (x.accession == "QCc:0000004"), runquality.qualityMetrics))), None) 
        metrics2add = list()
        if spectrum_acquisition_metrics_MS1:
            if starttime:
                metrics2add.extend(idfree.calcMSFrequency(spectrum_acquisition_metrics_MS1, start_time=starttime.value, ms_level=1))
                metrics2add.extend(idfree.calcMSDensity(spectrum_acquisition_metrics_MS1, start_time=starttime.value, ms_level=1))
            metrics2add.extend(idfree.calcSNMetrics(spectrum_acquisition_metrics_MS1, ms_level=1))
            metrics2add.extend(idfree.calcTrapTimeMetrics(spectrum_acquisition_metrics_MS1, ms_level=1))
            metrics2add.extend(idfree.calcMSSpread(spectrum_acquisition_metrics_MS1, ms_level=1))
        if spectrum_acquisition_metrics_MS2:
            if starttime:
                metrics2add.extend(idfree.calcMSFrequency(spectrum_acquisition_metrics_MS2, start_time=starttime.value, ms_level=2))
                metrics2add.extend(idfree.calcMSDensity(spectrum_acquisition_metrics_MS2, start_time=starttime.value, ms_level=2))
            metrics2add.extend(idfree.calcSNMetrics(spectrum_acquisition_metrics_MS2, ms_level=2))
            metrics2add.extend(idfree.calcTrapTimeMetrics(spectrum_acquisition_metrics_MS2, ms_level=2))
            metrics2add.extend(idfree.calcMSSpread(spectrum_acquisition_metrics_MS2, ms_level=2))
        if tandem_spectrum_metrics:
            metrics2add.extend(idfree.calcAnalysedSignalMetrics(tandem_spectrum_metrics, isolation_window_metrics))
            metrics2add.extend(idfree.calcPrecursorIntensityMetrics(tandem_spectrum_metrics))
            metrics2add.extend(idfree.calcPrecursorChargeMetrics(tandem_spectrum_metrics))
        if tic_tab:
            metrics2add.extend(idfree.calcESIStability(tic_tab))
            metrics2add.extend(idfree.calcTICSpread(tic_tab))
        runquality.qualityMetrics.extend(metrics2add)

    # m = [len(x.qualityMetrics) for x in MQR.run_qualities.values()]
    # logging.warning("{n} RunQualities with {mi} to {ma} metrics each.".format(n=len(MQR.run_qualities), mi=min(m), ma=max(m) ) )
    global rqs
    rqs = list(MQR.run_qualities.values())  # dict_values is a view , need to unpack

    # m = [len(x.qualityMetrics) for x in rqs]
    # logging.warning("{n} RunQualities with {mi} to {ma} metrics each.".format(n=len(rqs), mi=min(m), ma=max(m) ) )
    # for rq in rqs:
    #     for met in rq.qualityMetrics:
    #         if isinstance(met.value, {}.values().__class__):
    #             met.value = list(met.value)
    #             logging.warning("Converted '{}' metric values from dict_values".format(met.name))

    # import pickle
    # pickle.dump( MQR, open( "save.p", "wb" ) )
    # # mzqcfile = pickle.load( open( "save.p", "rb" ) )
    # with open('mqqcout.mzQC', 'w') as fh:
    #     fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs), readability=2))

    #TODO then set metrics according to the experimentalGrouping
    experiment_factors = set(MQR.experimentalGrouping.columns[MQR.experimentalGrouping.columns.str.startswith('Factor')])
    # TODO ask/let be chosen which factor but move functionality to get  grouped runs to MQanalysisResult
    considered_factor = next(iter(experiment_factors))
    experiment_groups = set(MQR.experimentalGrouping[considered_factor])
    grouped_runs = {group: set(MQR.experimentalGrouping['run name'][MQR.experimentalGrouping[considered_factor] == group]) for group in experiment_groups}
    
    # TODO collect set metric metric contributors
    idfm = ['Number of MS1 spectra', 
        'Number of MS2 spectra',
        'Fastest frequency for MS level 1 collection',
        'Fastest frequency for MS level 2 collection',
        'Slowest frequency for MS level 1 collection',
        'Slowest frequency for MS level 2 collection',
        'signal fall (10x) count',
        'signal jump (10x) count',
        'Median precursor charge in all MS2',
        'Charged spectra ratio +1 over +2',
        'Charged spectra ratio +3 over +2',
        'Charged spectra ratio +4 over +2']
    idbm = ['Number of identified spectra',
        'Number of identified peptides',
        'ID ratio',
        'Mean precursor charge in all MS2',
        'Explained base peak intensity median']
    # TODO idsplit (mostly from distros, evtl unavail. like S/N id vs unid)

    fm = pandas.DataFrame(columns=idfm+['factor','run']) 
    bm = pandas.DataFrame(columns=idfm+['factor','run'])
    for g,runs in grouped_runs.items():
        for run in runs:
            fmv = {metric_name: next(iter(list(filter(lambda x: (x.name == metric_name), MQR.run_qualities[run].qualityMetrics))), None) for metric_name in idfm}
            bmv = {metric_name: next(iter(list(filter(lambda x: (x.name == metric_name), MQR.run_qualities[run].qualityMetrics))), None) for metric_name in idbm}
            fmv.update({'factor': g, 'run': run})
            bmv.update({'factor': g, 'run': run})
            fm = fm.append(fmv, ignore_index=True)
            bm = bm.append(bmv, ignore_index=True)

    fm = fm.applymap(lambda x: x.value if isinstance(x,qc.QualityMetric) else x)
    bm = bm.applymap(lambda x: x.value if isinstance(x,qc.QualityMetric) else x)

    fg = fm.drop('run', 1).groupby('factor')
    between = (fg.count() * (fg.mean() - fm.mean())**2).sum() / (fm.run.size - len(grouped_runs))  # divided by number of observations - number of groups
    # ((f.drop(['run','factor'], 1) - f.drop('run', 1).groupby('factor').transform('mean'))**2).sum()  # f is the df with all None cols dropped too for quick exp.
    within = ((fm.drop(['run','factor'], 1) - fm.drop('run', 1).groupby('factor').transform('mean'))**2).sum() / (len(grouped_runs) -1) # sum of squares: sigma(value - group mean)**2 per group, all summed and divided by number of groups -1 
    
    # TODO once more for bm (and  then also for the new idsplit df)
    for g,runs in grouped_runs:
        mzqc_set = qc.SetQuality(qc.MetaDataParameters(label = g))
        # group variance
        
        # per group IsoFor
    # /sDMX outliers
    
    global_variables(output, zip)
    finale()
    # qccalculator maxquant --mqpar https://uk1s3.embassy.ebi.ac.uk/dda-mzqc-v2/test.PXD004364.mqpar.xml --zipuri https://uk1s3.embassy.ebi.ac.uk/dda-mzqc-v2/test.PXD004364.zip --eduri https://uk1s3.embassy.ebi.ac.uk/dda-mzqc-v2/test.PXD004364.sdrf.txt --mzmlfolder https://uk1s3.embassy.ebi.ac.uk/dda-mzqc-v2/ && echo $(date +%H:%M)

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Calculate quality metrics for given files.
    Multiple files input is possible (each after a "full/basic" COMMAND).
    All metrics of one QCCalculator execution will be stored in on output file.
    If you need separate mzQC files, please execute separately.
    For more information on the different COMMAND types, try qccalculator COMMAND --help
    """
    pass

cli.add_command(basic)
cli.add_command(full)
cli.add_command(maxquant)

if __name__ == "__main__":
    cli()
