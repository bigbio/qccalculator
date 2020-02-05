#!/usr/bin/env python
import os
import click
import logging
from datetime import datetime
import gzip
import pyopenms as oms
from mzqc import MZQCFile as qc

from .qccalculator import getBasicQuality, getIDQuality

rqs = list()
sqs = list()
out = str()
zp = False

#@click.pass_context
def finale():
    logging.warn("Calculated metrics from {} different input peak files".format(len(rqs)))
    logging.warn("Attempting to write results to {}{}".format(out, ".gz" if zp else ""))
    if any(rqs) or any(sqs) or out:
        if zp:
            with gzip.GzipFile(out+'.gz', 'w') as fh:
                fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs, out), readability=1).encode('utf-8'))
        else:
            with open(out, 'w') as fh:
                fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs, out), readability=0))
    logging.warn("Done. Thank you for choosing QCCalculator!")

def mzqc_assembly(rqs, sqs, out):
    # TODO check all the metrics to see which ontologies were used
    cv_qc = qc.ControlledVocabulary(ref="QC", 
        name="Proteomics Standards Initiative Quality Control Ontology",
        version="0.1.0", uri="https://github.com/HUPO-PSI/qcML-development/blob/master/cv/v0_1_0/qc-cv.obo")
    cv_ms = qc.ControlledVocabulary(ref="MS", 
        name="Proteomics Standards Initiative Mass Spectrometry Ontology",
        version="4.1.7", uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo")

    return qc.MzQcFile(version="0.1.0", 
                        creationDate=datetime.now().isoformat(), 
                        runQualities=rqs, 
                        setQualities=sqs, 
                        controlledVocabularies=[cv_qc, cv_ms])

@click.group(chain=True)
@click.option('--output', required=True, type=click.Path(), default="/tmp/out.mzQC", help="The path and name of the desired output file.")
@click.option('--zip/--no-zip', default=False, help="Apply gzip to the output. Appends '.gz' to the target filename and pretty formatting.")
def start(output, zip):
    """Calculate quality metrics for given files. 
       Multiple files input is possible (each after a "full/basic" COMMAND).
       All metrics of one QCCalculator execution will be stored in on output file. 
       If you need separate mzQC files, please execute separately.
       For more information on the different COMMAND types, try QCCalculator COMMAND --help"""
    logging.warn("Recieved output destination {}".format(output))
    global out 
    out = output
    if zip:
        global zp 
        zp = True

@start.command()
@click.argument('filename', type=click.Path(exists=True))
@click.option('--mzid', required=True, type=click.Path(exists=True), help="If you have a corresponding mzid file you need to pass it, too. Mutually exclusive to idxml.")
@click.option('--idxml', required=True, type=click.Path(exists=True), help="If you have a corresponding idxml file you need to pass it, too. Mutually exclusive to mzid.")
def full(filename, mzid=None): 
    """Calculate all possible metrics for these files. These data sources will be included in set metrics."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(click.format_filename(filename), exp)
    rq = getBasicQuality(exp)  
    
    ms2num = 0
    for x in rq.qualityMetrics:
        if x.name == "Number of MS2 spectra":
            ms2num = x.value

    if ms2num <1:
        logging.warn("We seem to have found no MS2 spectra which is unlikely to be true since you have also given some identifications. \
                We continue with symbolic value of 1 for the number of MS2 spectra, \
                however this means some metrics will invariably be incorrect!\
                Please make sure, we have the right inputs.")
        ms2num = 1

    if idxml and mzid:
        logging.warn("Sorry, you can only give one id file. Please choose one.")
        # TODO find a way to do this the click pallet way instead of warn
        return

    pros = list()
    peps = list() 
    if mzid:
        oms_id = oms.MzIdentMLFile()
        idf = mzid
    if idxml:
        oms_id = oms.MzIdentMLFile()
        idf = idxml
    if idf:
        oms_id.load(click.format_filename(idf), pros, peps)
        rq.qualityMetrics.extend(getIDQuality(exp, pros, peps, ms2num))
    rqs.append(rq)

    finale()
    
@start.command()
@click.argument('filename', type=click.Path(exists=True))
def basic(filename):
    """Calculate the basic metrics available from virtually every mzML file."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(click.format_filename(filename), exp)
    rq = getBasicQuality(exp)  
    rqs.append(rq)

    finale()

if __name__ == "__main__":
    start()
#QCCalculator --output cli-test.mzqc full --mzid tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.mzid tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.trfr.t3.mzML