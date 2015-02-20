Documentation for R analysis of pickle box microbes

James Meadow

February 19, 2015

------


This directory contains all R analysis for the second experiment of the pickle box manuscript. 

It all happens in several different scripts: 

* `cleanup.Rmd` = the first step. Input biom file, and then analyze/remove contamination. 
* `pb2014.Rmd` = this is where the map gets filled out, and data structured for analysis. 
* `finalAnalysis.Rmd` = where most of the actual analysis and figures are created. All information from the previous scripts is put together here. 


Several other parts that went into the manuscript: 

* `particle` = folder that contains analysis of particle count data from the experiment. this gets incorporated in the `analyzeSubsets.Rmd` script. 
* `quality` = a look at the quality from the 2 runs of paired reads. 
* `chamberConditions` = a quick visualization of the air flow rates and temp/relative humidity. These were within our target ranges, and are just used for quality checking. 
* `scripts` = a set of shell and perl scripts that were used to process the raw sequence data. The output of those scripts is a biom table that gets used by R for analysis. 


