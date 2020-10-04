# QCCalculator 
This folder contains classes and methods to calculate (_qccalculator_) QC metrics and create corresponding MZQC objects and their visualisation (_qcplots_).

## Calculation of metrics
The set of basic metrics are created by getBasicQuality which will consume at least one MS run and will encapsulate the created basic metrics in an MZQC.RunQuality. Other metric calculating methods however are producing MZQC.QualityMetrics and care must be taken to add these to the correct RunQuality elements.

### Metrics overview

### Aggregation metrics overview

#### basic
Spectrum_acquisition_metrics_MS1 (Spectrum acquisition metric values - MS1): RT, SN, peakcount, int
spectrum_acquisition_metrics_MS2 (Spectrum acquisition metric values - MS2): -""-
spectrum_topn (Spectra topn ranks): RT, rank
tandem_spectrum_metrics_MS2 (Tandem spectrum metric values - MS2): RT, precursor_int, surveyscan_int, precursor_err, precursor_mz, precursor_c, surveyscan_max
trap_metrics_MS1 (Trap metric values - MS1): RT, iontraptime
trap_metrics_MS2 (Trap metric values - MS2): RT, iontraptime
tic_tab (Total ion current chromatogram): RT, int
#### id
identification_accuracy_metrics (Identifications accuracy metric values): RT, MZ, delta_ppm, abs_err
identification_scoring_metrics (Identification scoring metric values): RT, c, score
identification_sequence_metrics (Identifications sequence metric values): RT, peptide, target
hydrophobicity_metrics (Hydrophobicity metric values): RT, gravy 
