#!/usr/bin/env python
from os.path import basename

import click
import logging
from datetime import datetime
import gzip
from typing import List

import pyopenms as oms
from click import command
from mzqc import MZQCFile as qc
from qccalculator import basicqc, idqc, idqcmq

rqs: List[qc.RunQuality] = list()
sqs: List[qc.SetQuality] = list()
out = str()
zp = False

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


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
  logging.warn("Received output destination {}".format(output))
  global out
  out = output
  if zip:
    global zp
    zp = True


def mzqc_assembly(rqs, sqs, out):
  """
  Setup the latest version of qcML ontology file

  Parameters
  ----------
  rqs: Run qualities
  sqs: Set Qualities
  out: output file file

  Returns
  -------

  """
  cv_qc = qc.ControlledVocabulary(ref="QC",
                                  name="Proteomics Standards Initiative Quality Control Ontology",
                                  version="0.1.0",
                                  uri="https://github.com/HUPO-PSI/qcML-development/blob/master/cv/v0_1_0/qc-cv.obo")
  cv_ms = qc.ControlledVocabulary(ref="MS",
                                  name="Proteomics Standards Initiative Mass Spectrometry Ontology",
                                  version="4.1.7", uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo")

  return qc.MzQcFile(version="0.1.0",
                     creationDate=datetime.now().isoformat(),
                     runQualities=rqs,
                     setQualities=sqs,
                     controlledVocabularies=[cv_qc, cv_ms])


# @click.pass_context
def finale():
  """
  Gunzip and compress the output files.

  Returns
  -------
  """
  logging.warn("Calculated metrics from {} different input peak files".format(len(rqs)))
  logging.warn("Attempting to write results to {}{}".format(out, ".gz" if zp else ""))
  if any(rqs) or any(sqs) or out:
    if zp:
      with gzip.GzipFile(out + '.gz', 'w') as fh:
        fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs, out), readability=1).encode('utf-8'))
    else:
      with open(out, 'w') as fh:
        fh.write(qc.JsonSerialisable.ToJson(mzqc_assembly(rqs, sqs, out), readability=0))
  logging.warn("Done. Thank you for choosing QCCalculator!")


def common_options(function):
  function = click.option('--mzml')(function)
  function = click.option('--output', required=True, type=click.Path(), default="out.mzQC",
                          help="The path and name of the desired output file.")(function)
  function = click.option('--zip', default=False,
                          help="Apply gzip to the output. Appends '.gz' to the target filename and pretty formatting.",
                          is_flag=True)(function)
  return function


@click.command('basic', short_help='Compute Qc metrics for an mzML file')
@common_options
def basic(mzml, output, zip):
  """
  Calculate the basic metrics available from virtually every mzML file.

  Parameters
  ----------
  mzml: compute the Qc metrics for a mzML file

  Returns
  -------

  """
  exp = oms.MSExperiment()
  oms.MzMLFile().load(click.format_filename(mzml), exp)
  rq = basicqc.getBasicQuality(exp)
  rqs.append(rq)

  global_variables(output, zip)
  finale()


@click.command("full", short_help='Compute the QC metrics for spectra file (mzML) and identification data')
@click.option('--mzid', type=click.Path(exists=True),
              help="If you have a corresponding mzid file you need to pass it, too. Mutually exclusive to idxml.")
@click.option('--idxml', type=click.Path(exists=True),
              help="If you have a corresponding idxml file you need to pass it, too. Mutually exclusive to mzid.")
@common_options
def full(mzid=None, idxml=None, mzml=None, output=None, zip=None):
  """
  Calculate all possible metrics for these files. These data sources will be included in set metrics.

  Parameters
  ----------
  mzid: the identification file in mzIdentML to compute the id qc metrics
  idxml: the identification file in idXML to compute the id qc metrics
  mzml: mzML to compute the qc metrics

  Returns
  -------

  """
  exp = oms.MSExperiment()
  oms.MzMLFile().load(click.format_filename(mzml), exp)
  rq = basicqc.getBasicQuality(exp)

  if idxml and mzid:
    with click.Context(command) as ctx:
      logging.warn("Sorry, you can only give one id file. Please choose one.")
      click.echo(command.get_help(ctx))
      return
  elif not idxml and not mzid:
    logging.warn("Sorry, you must give one id file in this mode.")
    click.echo(command.get_help(ctx))
    return

  ms2num = 0
  for x in rq.qualityMetrics:
    if x.name == "Number of MS2 spectra":
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
    rq.qualityMetrics.extend(idqc.getIDQuality(exp, pros, peps, ms2num))
  rqs.append(rq)

  global_variables(output, zip)
  finale()


@click.command("maxq", short_help="Compute the MaxQuant quality metrics")
@click.option('--zipurl', type=click.Path(exists=True), required=True,
              help="The URL to a max quant output zip file (must contain evidence.txt and parameters.txt).")
@click.option('--rawname', type=str, default="",
              help="The raw file name of interest (as in evidence.txt) without path or extension.")
@common_options
def maxq(zipurl, rawname, mzml, output, zip):
  """
  Calculate all possible metrics id and spectra for MaxQuant output. These data sources will be included in set metrics.

  Parameters
  ----------
  zipurl: zip file with all the MQ outputs.
  rawname: The files that will be use to compute the QC metrics.
  mzml: the mzML files to compute the QC metrics

  Returns
  -------

  """
  exp = oms.MSExperiment()
  oms.MzMLFile().load(click.format_filename(mzml), exp)
  rq = basicqc.getBasicQuality(exp)

  ms2num = 0
  for x in rq.qualityMetrics:
    if x.name == "Number of MS2 spectra":
      ms2num = x.value

  if ms2num < 1:
    logging.warn("We seem to have found no MS2 spectra which is unlikely to be true since you have also given some identifications. \
                We continue with symbolic value of 1 for the number of MS2 spectra, \
                however this means some metrics will invariably be incorrect!\
                Please make sure, we have the right inputs.")
    ms2num = 1

  try:
    mq, params = idqcmq.loadMQZippedResults(zipurl)
    if not rawname:
      logging.warning("Inferring rawname from mzML")
      rawname = basename(
        exp.getExperimentalSettings().getSourceFiles()[0].getNameOfFile().decode())  # TODO split extensions

    rq.qualityMetrics.extend(idqcmq.getMQMetrics(rawname, params, mq, ms2num))
    rqs.append(rq)
  except:
    logging.warn("Retrieving any results from the URL failed.")

  global_variables(output, zip)
  finale()


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
  """
  Calculate quality metrics for given files.
  Multiple files input is possible (each after a "full/basic" COMMAND).
  All metrics of one QCCalculator execution will be stored in on output file.
  If you need separate mzQC files, please execute separately.
  For more information on the different COMMAND types, try qccalculator COMMAND --help
  """


cli.add_command(basic)
cli.add_command(full)
cli.add_command(maxq)

if __name__ == "__main__":
  cli()
# qccalculator --output cli-test.mzqc full --mzid tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.mzid tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.trfr.t3.mzML
