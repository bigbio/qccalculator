#!/usr/bin/env python
from QCCalculator.qccalculator import getBasicQuality, getIDQuality
import os
import pyopenms as oms
from MZQC import MZQCFile as mzqc


def main():
    exp = oms.MSExperiment()
    oms.MzMLFile().load("tests/CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.trfr.mzML", exp)
    rq = getBasicQuality(exp)

if __name__ == "__main__":
    main()
