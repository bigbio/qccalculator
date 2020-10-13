import io
import zipfile
import urllib.request
import warnings
from itertools import chain
from typing import List, Tuple

import pandas
from Bio.SeqUtils import ProtParam
from mzqc import MZQCFile as mzqc

from qccalculator import utils

"""
Calculate id based metrics from MaxQuant result files
"""


def loadMQZippedResults(url: str) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
  """
    get_mq_zipped_evidence acquires the necessary MQ inputfiles from a URL to a zipped archive

    The predominant way identifications from MQ are stored in open mass spectrometry data repositories is in a zip file for the submission.
    The methods loads the archive and retrieves the QC metric relevant result files from the archive.

    Parameters
    ----------
    url : str
        A URL to a zip file with MQ result files.

    Returns
    -------
    Tuple[pandas.DataFrame,pandas.DataFrame]
        The parameters and evidence files from a MaxQuant result rehashed in accessible pandas dataframes.
    """
  with urllib.request.urlopen(url, timeout=10) as dl:
    with zipfile.ZipFile(io.BytesIO(dl.read())) as z:
      ef = 'evidence.txt'
      pf = 'parameters.txt'
      ld = dict()  # {'ev':'evidence.txt', 'pa':'parameters.txt'}
      dirs = dict()  # {f: {'ev':'evidence.txt', 'pa':'parameters.txt'} }
      for f in z.namelist():
        if z.getinfo(f).is_dir():
          dirs[f] = dict()
        elif f == ef:
          ld['ev'] = f
        elif f == pf:
          ld['pa'] = f

      if len(ld) < 2:
        for f in z.namelist():
          for d in dirs.keys():
            # exact expected match otherwise oddities like 'SEARCH/._parameters.txt' are picked up
            if f == d + ef:
              dirs[d]['ev'] = f
            elif f == d + pf:
              dirs[d]['pa'] = f

        dirs = {k: v for k, v in dirs.items() if len(v) > 0}
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
        ev = pandas.read_csv(e, sep='\t')
        ev.columns = map(str.lower, ev.columns)
      with z.open(ld['pa']) as p:
        pa = pandas.read_csv(p, sep='\t', dtype={'Parameter': str})
        pa.columns = map(str.lower, pa.columns)
        pa['parameter'] = pa['parameter'].str.lower()
        pa.set_index('parameter', inplace=True)
      return pa, ev


def getMQMetrics(target_raw: str, params: pandas.DataFrame, evidence: pandas.DataFrame, ms2num: int = 0) -> List[
  mzqc.QualityMetric]:
  """
    getMQMetrics calculates id based QC metrics from MaxQuant results as close as possible to the way they are calculated from regular id files.

    For a given raw file (name), the respective results are extracted from dataframes derived off the parameters and evidence files from a
    MaxQuant result (of potentially multiple raw files combined analysis). As many metrics similar or equal to those dependent of regular id files
    are calculated.

    Parameters
    ----------
    target_raw : str
        The name of the raw file (as per MaxQuant usage without file type extension)
    params : pandas.DataFrame
        Dataframe with data from the parameters result file as produced by MaxQuant and stratified column names
    evidence : pandas.DataFrame
        Dataframe with data from the evidence result file as produced by MaxQuant and stratified column names
    ms2num : int, optional
        The total number of tandem spectra as from the id-free metrics, by default 0

    Returns
    -------
    List[mzqc.QualityMetric]
        A list of QualityMetrics close to what is calculated from a regular id-based QC calculation.
    """
  if not target_raw in evidence['raw file'].unique():
    return list()  # TODO warn
  else:
    mq_metrics: List[mzqc.QualityMetric] = list()
    # https://stackoverflow.com/questions/17071871/how-to-select-rows-from-a-dataframe-based-on-column-values
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
    psims = utils.obtainOntology("psi-ms")

    name_indexed = {psims[x].name: psims[x] for x in psims}
    score_indexed = {x.name: x for x in chain(psims['MS:1001143'].subclasses(), psims['MS:1001153'].subclasses(),
                                              psims['MS:1002347'].subclasses(), psims['MS:1002363'].subclasses())}

    if score_type in name_indexed:
      if not score_type in score_indexed:
        warnings.warn("Score type does not correspond to a score type in the OBO, proceed at own risk.", Warning)
        score_col_name = name_indexed[score_type].id
      else:
        score_col_name = score_indexed[score_type].id
    else:
      warnings.warn("OBO does not contain any entry matching the identification score, proceed at own risk.", Warning)
      score_col_name = score_type

    identification_scoring_metrics = target_mq[['retention time', 'charge', 'score']].rename(
      columns={'retention time': 'RT', 'charge': 'c', 'score': score_type}).to_dict(orient='list')
    mq_metrics.append(
      mzqc.QualityMetric(cvRef="QC",
                         accession="QC:0000000",
                         name="Identification scoring metric values",
                         value=identification_scoring_metrics)
    )

    # TODO comparison column with qccalculator dppm values
    # TODO RT/native id?
    identification_accuracy_metrics = target_mq[['ms/ms m/z', 'mass error [ppm]', 'uncalibrated mass error [da]']] \
      .rename(columns={'ms/ms m/z': 'MZ', 'mass error [ppm]': 'delta_ppm', 'uncalibrated mass error [da]': 'abs_error'})
    identification_accuracy_metrics['abs_error'] = identification_accuracy_metrics['abs_error'].abs()
    identification_accuracy_metrics = identification_accuracy_metrics.to_dict(orient='list')
    mq_metrics.append(
      mzqc.QualityMetric(cvRef="QC",
                         accession="QC:0000000",
                         name="Identifications accuracy metric values",
                         value=identification_accuracy_metrics)
    )

    hydrophobicity_metrics = target_mq[['retention time', 'sequence']].rename(
      columns={'retention time': 'RT', 'sequence': 'peptide'})
    hydrophobicity_metrics['gravy'] = hydrophobicity_metrics['peptide'].apply(
      lambda x: ProtParam.ProteinAnalysis(x).gravy())
    hydrophobicity_metrics = hydrophobicity_metrics[['RT', 'gravy']].to_dict(orient='list')
    mq_metrics.append(
      mzqc.QualityMetric(cvRef="QC",
                         accession="QC:0000000",
                         name="Hydrophobicity metric values",
                         value=hydrophobicity_metrics)
    )

    # TODO target/decoy info available??
    identification_sequence_metrics = target_mq[['sequence', 'retention time', 'ms/ms scan number']].rename(
      columns={'sequence': 'peptide', 'retention time': 'RT', 'ms/ms scan number': 'native_id'}).to_dict(orient='list')
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
                         value=float(len(target_mq)) / float(ms2num))
    )

    return mq_metrics
