import re
from collections import defaultdict
from typing import Any,  Dict, List, Tuple, Pattern

from Bio import  SeqRecord
from mzqc import MZQCFile as mzqc
import pyopenms as oms
import numpy as np

from qccalculator import utils

"""
Methods to calculate quality metrics related to digestion enzyme use during sample preparation and identification process
"""


def getCoverageRatios(pro_ids: oms.ProteinIdentification,
                      pep_ids: List[oms.PeptideIdentification],
                      fasta=Dict[str, SeqRecord.SeqRecord], fetch=False) -> List[mzqc.QualityMetric]:
  """
    getCoverageRatios calculates the coverage ratios per protein from the identified searchspace.

    Calculating the coverage from all individual petide identification also requires all protein
    sequences expected to be known. For this there are two options, either retrieve the sequences
    from the originally used fasta, or try to retrieve the sequences via UniProt through the
    accessions with the PeptideHits.

    Parameters
    ----------
    pro_ids : List[oms.ProteinIdentification]
        The PyOpenMS ProteinIdentification as from reading a common identification file
    pep_ids : List[oms.PeptideIdentification]
        List of PyOpenMS PeptideIdentification as from reading a common identification file
    fasta : [type], optional
        Dictionary of sequences from a fasta file, Dict[accession,SeqRecord] by default Dict[str,SeqRecord.SeqRecord]
    fetch : bool, optional
        If set true, will attempt to retrieve sequences by accession, is ignored if `fasta` is provided, by default False

    Returns
    -------
    List[mzqc.QualityMetric]
        [description]
    """
  metrics: List[mzqc.QualityMetric] = list()

  # check all proteinhits have seq set
  # first via proteinhits, missing then either via fasta or
  # calc coverage
  missing_acc = list()
  nup = list()
  for p in pro_ids.getHits():
    ac = p.getAccession()
    nup.append(oms.ProteinHit(p))
    if not p.getSequence():
      if fasta:
        nup[-1].setSequence(str(fasta.get(ac, SeqRecord.SeqRecord('')).seq))
      # if still no sequence
      if not p.getSequence():
        missing_acc.append(ac)

  if missing_acc:
    uniprot = {x.id: x for x in utils.getUniProtSequences(missing_acc)}
    for n in nup:
      ac = n.getAccession()
      if not n.getSequence():
        n.setSequence(str(uniprot.get(ac, SeqRecord.SeqRecord('')).seq))
    urx = re.compile('\w*\|(\w*)\|\w*')
    uniprot = {re.search(urx, x.id).group(): x for x in utils.getUniProtSequences(missing_acc)}
    del uniprot['']
    for n in nup:
      ac = n.getAccession()
      if not n.getSequence():
        n.setSequence(str(uniprot.get(ac, SeqRecord.SeqRecord('')).seq))

  coverage_tab: Dict[str, List[Any]] = defaultdict(list)
  na = [n.getAccession() for n in nup if not n.getSequence()]
  nup = [n for n in nup if n.getSequence()]

  pro_ids.setHits(nup)
  pro_ids.computeCoverage(pep_ids)

  for p in pro_ids.getHits():
    coverage_tab['Accession'].append(p.getAccession())
    coverage_tab['Coverage'].append(p.getCoverage())
    coverage_tab['Length'].append(len(p.getSequence()))
    # TODO figure out decoy string by fasta
    coverage_tab['TD'].append('decoy' if 'decoy' in p.getAccession().lower() else 'target')

  for n in na:
    coverage_tab['Accession'].append(n.getAccession())
    coverage_tab['Coverage'].append('NA')
    coverage_tab['Length'].append('NA')
    coverage_tab['TD'].append('decoy' if 'decoy' in n.getAccession().lower() else 'target')

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Protein coverage",
                       value=coverage_tab)
  )

  return metrics


def getPeptideLengthMetrics(identification_sequence_metrics: mzqc.QualityMetric) -> List[mzqc.QualityMetric]:
  """
    describePeptideLengthMetrics calculates the descriptive statistics metrics for identified sequences' length

    From the proto-metrics on identification sequences, the function calculates descriptive statistics metrics for
    the distribution of peak density from all involved mass spectra.
    Namely, mean, standard deviation, Quartiles, and 1.5*IQR outliers.

    Parameters
    ----------
    identification_sequence_metrics : mzqc.QualityMetric
        QualityMetric with 'peptide' value, filtered for final outcome

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  regex_mod = r'(\([^\(]*\))'
  regex_noaa = r'([^A-Za-z])'
  # TODO test this: '.(iTRAQ4plex)M(Oxidation)C(Carbamidomethyl)HNVNR'
  lengths = np.array(
    [len(re.sub(regex_noaa, '', re.sub(regex_mod, '', x))) for x in identification_sequence_metrics.value['peptide']])

  q1, q2, q3, s, m, ol = utils.extractDistributionStats(lengths)
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Identified peptide lengths Q1, Q2, Q3",
                                    value=[q1, q2, q3])
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Identified peptide lengths sigma",
                                    value=s)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Identified peptide lengths mean",
                                    value=m)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Identified peptide lengths +/-1.5*IQR outlier",
                                    value=ol)
                 )

  return metrics


def matchEnzyme(enzre: Pattern, pepseq: str) -> Tuple[int, int]:
  """
    matchEnzyme matches a peptide sequence against a regular expression pattern

    The regular expression pattern is matched against the peptide sequence
    and the matches counted according to their position/type. For that, the
    peptide sequences should include the amino acids before and after as they
    occurr in the protein. It is distinguished between match/semi (at the peptide
    end(s)) and internal matches. This makes for the return values of
    match/semi/none and number of internal matches (a.k.a. missed cleavages).

    Parameters
    ----------
    enzre : re._pattern_type
        Pattern as created by re.compile(...)

    pepseq : str
        amino acid sequence string

    Returns
    -------
    Tuple(int,int)
        First int indicates matched (2) semi matched (1) none (0), second int is the count or internal matches
    """
  matches = np.array([x.start() if x.start() == x.end() else None for x in enzre.finditer(pepseq)])
  is_matched = False
  is_semi = False
  if 1 in matches and len(pepseq) - 1 in matches:
    is_matched = True
    internal_matches = len(matches) - 2
  elif 1 in matches or len(pepseq) - 1 in matches:
    internal_matches = len(matches) - 1
    is_semi = True
  else:
    internal_matches = len(matches)

  return 2 if is_matched else 1 if is_semi else 0, internal_matches


def getEnzymeContaminationMetrics(pep, pro, force_enzymes=False) -> List[mzqc.QualityMetric]:
  """
    getEnzymeContaminationMetrics calculates enzyme and enzyme contamination metrics from the
    identifications given.

    The function calculates the number of missed cleavages (internal), peptide length distribution,
    and peptide boundaries matching known enzyme patterns from the given identifications. Matching
    against digestion enzyme patterns other than the enyme used for identification processess has to
    be switched with 'force_enzymes' and is sensible if the identification was conducted with
    unspecific cleavage to detect enzyme contamination or enzyme setting mixup is suspected.

    Parameters
    ----------
    pro : List[oms.ProteinIdentification]
        List of PyOpenMS ProteinIdentification as from reading a common identification file
    pep : List[oms.PeptideIdentification]
        List of PyOpenMS PeptideIdentification as from reading a common identification file
    force_enzymes : bool, optional
        If set, will force checking the identified peptide sequences against other known
        digestion enzyme patterns. By default False

    Returns
    -------
    List[mzqc.QualityMetric]
        List of resulting QualityMetrics
    """
  metrics: List[mzqc.QualityMetric] = list()

  # include all psm actually does not make much sense to assess the enzyme efficiency
  gre = {pro[0].getSearchParameters().digestion_enzyme.getName():
           re.compile(pro[0].getSearchParameters().digestion_enzyme.getRegEx())}

  # TODO pyopenms wrappers for DigestionEnzymeDB etc
  # li: List = list()
  # oms.DigestionEnzymeDB().getAllNames(li)
  # ore = {e: re.compile(oms.DigestionEnzymeDB().getEnzyme(e).getRegEx()) for e in li
  #            if e not in gre and e != 'no cleavage'}

  enzymematch_tab: Dict[str, List[Any]] = defaultdict(list)
  missed_ranks = list()
  matched_ranks = list()
  # alt = dict()
  for i, pepi in enumerate(pep):
    pepi.sort()
    spec_id = pepi.getMetaValue('spectrum_reference') \
      if pepi.metaValueExists('spectrum_reference') else i
    for i, h in enumerate(pepi.getHits()):
      pepseq = h.getPeptideEvidences()[0].getAABefore() \
               + h.getSequence().toUnmodifiedString() \
               + h.getPeptideEvidences()[0].getAAAfter()

      is_matched, internal_matches = matchEnzyme(next(iter(gre.values())), pepseq)
      if i == 0:
        enzymematch_tab['native_id'].append(spec_id)
        enzymematch_tab['matched'].append(is_matched)
        enzymematch_tab['missed'].append(internal_matches)
      else:
        missed_ranks.append(internal_matches)
        matched_ranks.append(is_matched)

      # if force_enzymes or not is_matched:
      #     oth_enz_matched = {k: matchEnzyme(v, pepseq) for k,v in ore.items()}
      #     alt[spec_id] = oth_enz_matched

  if len(missed_ranks):
    arr = np.array(missed_ranks)
    q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)
    metrics.append(mzqc.QualityMetric(cvRef="QC",
                                      accession="QC:0000000",
                                      name="Q1, Q2, Q3 of missed clevage counts for all lower rank identifications.",
                                      value=[q1, q2, q3])
                   )

    metrics.append(mzqc.QualityMetric(cvRef="QC",
                                      accession="QC:0000000",
                                      name="Sigma of missed clevage counts for all lower rank identifications.",
                                      value=s)
                   )

    metrics.append(mzqc.QualityMetric(cvRef="QC",
                                      accession="QC:0000000",
                                      name="Mean of missed clevage counts for all lower rank identifications.",
                                      value=m)
                   )

    metrics.append(mzqc.QualityMetric(cvRef="QC",
                                      accession="QC:0000000",
                                      name="Missed clevage count for all lower rank identifications +/-1.5*IQR outlier",
                                      value=ol)
                   )

  if len(matched_ranks):
    mdl: Dict[int, int] = defaultdict(int)
    arr = np.array(matched_ranks)
    uniq, counts = np.unique(arr, return_counts=True)
    mdl.update(dict(zip(uniq, counts)))
    metrics.append(
      mzqc.QualityMetric(cvRef="QC",
                         accession="QC:0000000",
                         name="Match/semi/none counts for all lower rank identifications.",
                         value=[mdl[2], mdl[1], mdl[0]])
    )

  metrics.append(
    mzqc.QualityMetric(cvRef="QC",
                       accession="QC:0000000",
                       name="Missed cleavages",
                       value=enzymematch_tab)
  )

  arr = np.array(enzymematch_tab['missed'])
  q1, q2, q3, s, m, ol = utils.extractDistributionStats(arr)
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Q1, Q2, Q3 of missed clevage counts for top identifications.",
                                    value=[q1, q2, q3])
                 )
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Sigma of missed clevage counts for top identifications.",
                                    value=s)
                 )
  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Mean of missed clevage counts for top identifications.",
                                    value=m)
                 )

  metrics.append(mzqc.QualityMetric(cvRef="QC",
                                    accession="QC:0000000",
                                    name="Missed clevage count for top identifications +/-1.5*IQR outlier",
                                    value=ol)
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
# TODO unit tests
# XMAGHHHEHEQERDHEQEHEHDSLQRP
# KPNPASMX
# RPSTNDPTSCCSX
