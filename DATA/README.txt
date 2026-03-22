Input data
----------
This folder holds CSV files exported from Thermo TraceFinder, the peak
detection and compound identification software bundled with Thermo Fisher
Scientific GC-MS instruments. TraceFinder performs deconvolution and library matching against
NIST mainlib / replib, then exports one CSV per sample run.

This folder is excluded from version control (.gitignore).
Add your own CSV exports here before running the pipeline.


How to export from TraceFinder
------------------------------
Export one file per sample (including blanks). The resulting CSV starts
directly with the column-header row -- no additional preamble is required,
but the pipeline also handles a single non-data header line if present.


File naming convention
----------------------
Files are classified by their filename prefix as configured in config.py:

  Blank*.csv      -- process blanks  (e.g. Blank1_run1.csv)
  <prefix>*.csv   -- biological samples, one group per prefix

Example (default config with two sample groups):
  Blank1_2026_3_4.csv
  S1_2026_3_4.csv      -- group "S"   (prefix "S")
  S-R1_2026_3_4.csv    -- group "S-R" (prefix "S-R", matched first)

Put more specific prefixes before less specific ones in SAMPLE_GROUPS so
that "S-R" is matched before "S".


Required columns
----------------
The pipeline reads three columns from each CSV. Their names must match
exactly (TraceFinder uses these names by default):

  Retention Time   -- peak apex RT in minutes
  Reference m/z    -- accurate monoisotopic m/z from the Orbitrap detector
  Area             -- integrated peak area (use HEIGHT if preferred;
                      set VALUE_COL = "Height" in config.py)

All other columns exported by TraceFinder are ignored.


Example CSV (abbreviated)
-------------------------
TraceFinder exports a comma-separated file with one row per detected peak.
Unidentified peaks appear as "Peak@<RT>" in the Component Name column and
are retained by the pipeline just like named compounds.

Component Name,Retention Time,Reference m/z,Area,...,Sample Name,File name
"Chloroiodomethane",3.086,175.888397,60495308,...,S1_2026_3_4,C:\...\S1.raw
"3-Hexenal",4.044,41.038414,189949422,...,S1_2026_3_4,C:\...\S1.raw
"4-Penten-1-ol, 3-methyl-, acetate",10.621,67.054230,237184140,...,S1_2026_3_4,C:\...\S1.raw
"Indole",19.072,117.057228,916555386,...,S1_2026_3_4,C:\...\S1.raw
"Peak@21.166",21.166,87.044098,1095820,...,S1_2026_3_4,C:\...\S1.raw

The full header exported by TraceFinder contains these columns:
  Component Name, Retention Time, Reference m/z, Area, Height, TIC,
  Formula (mol ion), CAS No., SI, RSI, HRF Score, RHRF Score, Total Score,
  Selected Column Type, Calculated RI, Library RI, Delta RI,
  Library Name, Library ID Number, Sample Name, File name


Notes
-----
- RT units must be minutes (TraceFinder default).
- The "Reference m/z" column contains the accurate mass measured by the
  Orbitrap; this is used for optional m/z-based feature sub-clustering
  (USE_MZ = True) and for RT alignment (MZ_ALIGN_TOLERANCE).
- Library matching quality columns (SI, RSI, scores) are not used by the
  pipeline -- all detected peaks are carried forward regardless of score.
- If a compound is detected more than once at the same RT in a single
  sample, the pipeline keeps the entry with the highest area.
