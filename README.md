## STABILITYCALC

### by Blaise Frederick and Daniel M. Drucker

Stabilitycalc is part of a suite of python tools we use to monitor the
stability of our fMRI system.  We use a version of the fBIRN procedure
to run daily phantom quality assurance tests on our research MR system.
These phantom scans are then converted to NIFTI format and analyzed 
with stabilitycalc to look for system instabilities that could degrade
our fMRI data.

The stabilitycalc package generates reports for individual scans, and
tabulates performance over time, generating summary reports that put
current performance into the context of normal operation.  We use this
to quickly identify scanner problems, and the individual reports give
a good indication of what sort of problem we are having, which in turn
helps us troubleshoot much more effectively.

## Requirements

Stabilitycalc checks for some of its requirements before running. In general,
it requires a standard numpy/scipy/matplotlib environment, as well as
`nibabel`, `nipype`, `dicom`, `pandas`, and `seaborn`.
