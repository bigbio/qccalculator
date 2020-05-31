import hashlib
import re
import requests
import warnings
from io import StringIO
import urllib
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

import pyopenms as oms
import numpy as np
import pronto
from Bio import SeqIO, SeqRecord

"""
Utility functions that do not contribute directly to QC calculations 
"""

def sha256fromfile(abs_file_path: str) -> str:
    """
    sha256fromfile will create a sha256 digest from the file at given path.

    To preserve memory and speed up the digest, 
    the file is digested with the help of a memoryview and hashlib.sha256().update.

    :raises 
    
    Parameters
    ----------
    abs_file_path : str
        The absolute path to the file to digest

    Returns
    -------
    str
        The cast or unchanged argument
    
    Raises
    ------
    FileNotFoundError 
        If abs_file_path is not a file
    """
    sha  = hashlib.sha256()
    b  = bytearray(128*1024)
    mv = memoryview(b)
    
    with open(abs_file_path, 'rb', buffering=0) as f:
        for n in iter(lambda : f.readinto(mv), 0):
            sha.update(mv[:n])
    return sha.hexdigest()

def cast_if_int(pot_int: Any) -> Union[int,Any]:
    """
    cast_if_int convenience function to cast to int
    
    Due to the frequent use of numpy.dtypes and pyOpenMS return of binary encode strings, 
    this function will ease the level of verbosity.

    Parameters
    ----------
    pot_int : Any
        The potential int value

    Returns
    -------
    Union[int,Any]
        In case the argument is cast-able into int, will return that int, unchanged argument otherwise.
    """
    try:
        return int(pot_int)
    except ValueError as e:
        return pot_int

def spec_native_id(spec: oms.MSSpectrum) -> Union[int,None]:
    """
    spec_native_id convenience function to retrieve the native id number from a spectrum

    Since the spectrums native id is a string formatted with much additional, albeit 
    usually redundant information, this method cuts through the clutter and extracts 
    the numerical id.

    Parameters
    ----------
    spec : oms.MSSpectrum
        Spectrum to get the native id from

    Returns
    -------
    Union[int,None]
        Return is None if spectrum native id cannot be interpreted (e.g. not of scan=number format)
    """
    spre = spec.getNativeID()
    if spre:
        matches = re.findall("scan=(\d+)$", spre)
        if len(matches)!=1:  # should really never be >1 with the `$`
            return None
        else:
            return cast_if_int(matches[0])
    else:
        return None

def pep_native_id(p: oms.Peptide) -> Union[int,None]:
    """
    pep_native_id convenience function to retrieve the native id number from an identification
    
    Counterpart to spec_native_id.
    Identifications loaded from mzid et al. should carry the native id to which spectra they 
    carry the identification information (as 'spectrum_reference'). Since the spectrums 
    native id is a string formatted with much additional, albeit usually redundant 
    information, this method cuts through the clutter and extracts the numerical id.

    Parameters
    ----------
    p : oms.Peptide
        PeptideIdentification from which to get the native id of the involved spectrum

    Returns
    -------
    Union[int,None]
        Return is None if native id cannot be interpreted (e.g. not of scan=number format)
    """
    spre = p.getMetaValue('spectrum_reference')
    if spre:
        matches = re.findall("scan=(\d+)$", spre)
        if len(matches)!=1:  # should really never be >1 with the `$`
            return None
        else:
            return cast_if_int(matches[0])
    else:
        return None

def getMassDifference(theo_mz: float, exp_mz: float, use_ppm: bool=True)-> float:
    """
    getMassDifference convenience function to easily switch the delta mass to either [ppm] or [Da] format.

    Given two masses, the calculated result will be the delta mass, in [ppm] if requested. 
    The difference is **not** absolute.

    Parameters
    ----------
    theo_mz : float
        First mass
    exp_mz : float
        Second mass
    use_ppm : bool, optional
        switch from simple [Da] difference to [ppm], by default True

    Returns
    -------
    float
        [description]
    """
    error: float = (exp_mz - theo_mz)
    if use_ppm:
        error = error / (theo_mz * 1e-6)
    return error

def getTrapTime(spec: oms.MSSpectrum, acqusition_unavailable=False) -> float:
    """
    getTrapTime for a given MSn spectrum, return the ion trap collection time applied during acquisition.

    The ion collection time, usually the same for spectra of one level from one
    run is taken from the mzML file. The value is sourced from 'MS:1000927' but 
    returned as a negative value if no cvTerm was available (not mandatory in mzML).
    This means in particular that MS1 spectra will all have negative values returned.
    In case pyopenms < 2.5.0 , use the acqusition_unavailable flag 
    and execute TrapTimeTool before QCCalculator execution.    

    Parameters
    ----------
    spec : oms.MSSpectrum
        Spectrum to get the trap time from 
    acqusition_unavailable : bool, optional
        In case access to AcquisitionInfo through pyopenms is unavailable, by default False

    Returns
    -------
    float
        Ion trap collection time in [ms] for given MSn
    """    
    tt = -1.0
    if acqusition_unavailable:
        if spec.metaValueExists('MS:1000927'):
            tt = spec.getMetaValue('MS:1000927')        
    else:
        if not spec.getAcquisitionInfo():
            for j in spec.getAcquisitionInfo():
                if j.metaValueExists("MS:1000927"):
                    tt = j.getMetaValue("MS:1000927")
                    break
    return tt

def extractDistributionStats(value_array: np.array) -> Tuple:
    """
    extractDistributionStats pulls descriptive distribution stats from an numpy array of values

    Extracted are the quartiles, sigma, mean, and outlier values ><1.5*IQR (in no guaranteed order)

    Parameters
    ----------
    value_array : np.array
        numpy array of ndim=1, all values of one type

    Returns
    -------
    Tuple
        In order the values Q1, Q2, Q3, sigma, mean, outliers
    """    
    q1, q2, q3 = np.quantile(value_array, [.25,.5,.75])
    s=np.std(value_array)
    m=np.mean(value_array)
    
    low_out = q1-(1.5*(q3-q1))
    high_out = q3+(1.5*(q3-q1))
    ol=np.extract((value_array<low_out) | (value_array>high_out), value_array)
    return q1, q2, q3, s, m, ol

def getUniProtSequences(accessions: List[str]) -> Union[List[SeqRecord.SeqRecord],None]:
    """
    getUniProtSequences retrieves the 

    The function will formulate a query to uniprot and parse the result into a list of Bio.SeqRecord s.
    
    No check on the completeness of the result is done here. 
    --------------------------------------------------------
    However, no isoforms are reported. So the number of SeqRecord s should be equal to the number of 
    accessions queried initially.

    Parameters
    ----------
    accessions : List[str]
        The list of uniprot accessions to query for their sequence

    Returns
    -------
    Union[List[SeqRecord.SeqRecord],None]
        The resulting list of SeqRecord accessions
    
    # https://docs.python.org/3/library/xml.etree.elementtree.html#pull-api-for-non-blocking-parsing
    """    
    acc = '+OR+'.join(['id:'+a for a in accessions])
    params = {"query": acc, "format": "fasta", "include": "no"}  # no isoforms
    response = requests.get("https://www.uniprot.org/uniprot/", params=params, verify=False)  # ugh, https certificate verification does not work OOTB with uniprot.org
    if not response.ok:
        response.raise_for_status()
        warnings.warn("UniProt query for sequences unsuccessful.")
        return None
    seqs = [s for s in SeqIO.parse(StringIO(response.text),format='fasta')]
    
    return seqs

def obtainOntology(url: str) -> pronto.Ontology:
    """
    obtainOntology provides pronto ontology objects and handles aspects of provision 

    An ontology can be requested by URL or common names which draws the latest available 
    version of the respective ontology. Available are 'psi-ms', 'psi-qc', 'units', 'pride', 'MS', 'QC', 'UO'.

    Parameters
    ----------
    url : str
        Either the url or a common name

    Returns
    -------
    pronto.Ontology
        pronto ontology
            
    Raises
    ------
    poof
        general NameError with url information
    """
    usual_suspects = {"psi-ms": "https://github.com/HUPO-PSI/psi-ms-CV/releases/download/4.1.31/psi-ms.obo",
                    "psi-qc": "https://raw.githubusercontent.com/HUPO-PSI/mzQC/master/cv/qc-cv.obo",
                    "units": "http://purl.obolibrary.org/obo/uo.owl",
                    "pride": "http://purl.obolibrary.org/obo/pride_cv.obo"}
    usual_suspects.update({"MS": usual_suspects["psi-ms"], 
                        "QC": usual_suspects["psi-qc"],
                        "UO": usual_suspects["units"]})
    # TODO move hardcoded urls to file
 
    if url in usual_suspects:
        url = usual_suspects[url]
    try:
        with urllib.request.urlopen(url, timeout=10) as obo_in:
            obo = pronto.Ontology(obo_in)
    except Exception as e:
        poof = NameError("Unable to obtain {}; please make sure the url exists, is available, and contains a parseable ontology.".format(url))
        raise poof from e
        # TODO separate 404 from connection/transmission error
    return obo
