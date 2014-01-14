This directory contains an analysis record for particle counts and 16S
sequence data.

The purpose of the experiment was to determine the nature and
identifiability of the human microbial cloud.

Within the `manuscript_code` folder:

* `contamination`: control samples were assessed for lab contamination,
and a few abundant contaminants were removed from further analysis.
* `cleanup`: put sequence data into shape for analysis.
* `analysis`: most ecological statistical anaysis for the manuscript.
* `taxon_analysis`: statistical steps related to individual bacterial
taxa - indicator analysis takes a bit of time, so it was put into its
own script.
* `tables_figures`: all tables, figures, and reported statistics come
from this script.
* `particles`: particle count data was analyzed separately, all within
this script.
* `data`: ~~~raw data as well as~~~ R output from individual scripts lives
here. Some raw data had to move to figshare (
http://dx.doi.org/10.6084/m9.figshare.900389)
 because the main OTU table was way too big. 
So download `otu_table_r.txt.gz` from there and unzip it into the `data` folder. 
