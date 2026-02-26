PyAstro Tools — Variable Star Observation Suite

A collection of Python desktop applications for variable star observers who submit data to the AAVSO (American Association of Variable Star Observers). All tools share a consistent GUI built with Tkinter and follow the AAVSO Extended File Format standard.


PyAAVSOGenerator

Converts raw photometry output from AstroImageJ (AIJ) into properly formatted AAVSO Extended Format reports ready for submission. Supports ensemble and single-comparison-star reductions, check-star validation, and generates light-curve plots alongside the report. Configurable per session for star name, filter, observer code, chart ID, and transformation flags.


PyAAVSOPlanner

An observation planning tool that computes and plots the altitude of a target star throughout a night for any given site and date. Fetches recent AAVSO observations of the target from the AAVSO web API and overlays them on the visibility chart, helping you decide when and whether to observe. Supports configurable site coordinates, altitude thresholds, and time-step sampling.


PyPeriodAnalysis

Performs Lomb-Scargle period analysis on AAVSO Extended Format data files. Load your photometry, set period search bounds, optionally filter by passband, and the tool computes the periodogram, identifies the strongest period candidate, and displays phase-folded light curves — all in an interactive GUI.


PyTransformApplier

Applies photometric transformation coefficients (BV or BVR mode) to AAVSO Extended Format files using coefficients stored in a VPhot.ini file. Outputs transformed data as .xlsx and .txt files in AAVSO Extended Format, and generates light-curve .png images — making it straightforward to produce publication-ready, color-corrected photometry.


PyVarstarMinimum

A minima calculator for eclipsing binary and pulsating variable stars. Looks up ephemeris data (epoch and period) directly from the VSX (Variable Star Index) database, then propagates the epoch forward or backward to list predicted minima/maxima closest to a user-specified date and time. Useful for scheduling eclipse observations or confirming O-C timing.


Requirements

Python 3.9+, astropy, pandas, numpy, matplotlib, requests (PyAAVSOPlanner), tkinter (usually bundled with Python on Windows/macOS).

Install dependencies with: pip install astropy pandas numpy matplotlib requests openpyxl


Author

Nikola Antonov
nikola.antonov@iaps.institute
https://astro.iaps.institute
