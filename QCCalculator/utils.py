import hashlib
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union

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
