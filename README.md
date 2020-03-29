# TreeSwallow_SURF

In this repository, there are two core R files, and several data files. The core R files include:

* `TreeSwallow_Rscript.R` - this is the bulk of the analysis, and is organized by paper section. 

* `CF_Tree_Swallows_Git.R` - this is a supplemental analysis done using causal forests. Results from this do not appear in the main text. 

The data files included are loaded at various points using my own computer's local directory. To reproduce results, these will have to be changed to whatever directory the data are downloaded into. Many of these files are .rda files, including processed results, but .shp GIS files are also included as a .zip file, for space purposes. These will need to be extracted before running the code.
