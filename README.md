# fundicalc
fundicalc is an R implementation that finds optimal Rho and branch length values given a FunDi split point in a tree.

Authors: Justin Daneau, Edward Susko, Hector Baños, Andrew Roger

Current Version: 0.5.0 (2022-09-15)

### Introduction
[todo]

### Installation
This software has only been tested to work on Linux and MacOS systems; Windows systems are not supported due to incompatibility. Ensure that R is installed on the system.

fundicalc comes bundled with two main R files (rtsplit,fundicalc), two supplemental R script packages (convert_tree.R, root_tree.R), and three binaries (rert, treecns, del-taxa-utreec). The binaries were written in the C language, originally for use with gfmix. 

In order to use the `fundicalc` script, there must be a valid installation of IQ-TREE on the system which can be called from any location. Currently, the version called in the script is `iqtree2203`. To change this, edit this line at the top of the script: `iqtree <- "iqtree2203"`

All files must be placed in the same directory for the scripts to function. Run the command `chmod a+x [filename]` for all of `fundicalc`, `del-taxa-utreec`, `rert`, `rtsplit`, and `treecns` to allow command-line execution of the files. 

### Program Instructions
There are two main files, `rtsplit` and `fundicalc`, which are responsible, respectively, for splitting up a tree into trees A, B and AB according to the FunDi split, and for running IQ-TREE to optimize the values of Rho and `t` (the split edge length) according to the formula `L_overall = (rho)(L_ab) + (1 - rho)(l_a)(l_b)`.

There are four required files required to run the software: a sequence file, a newick tree file, a "root" file indicating the set of taxa on one side of the root for the tree indicated in the tree file, and a "split" file indicating the set of taxa on one side of the FunDi split. The latter two files are of the same format, with space-separated integers on one line. 

To properly run fundicalc, rtsplit must first be run. This is because fundicalc relies on output from rtsplit. 

#### Input commands

- rtsplit
The command for rtsplit looks like `rtsplit -s seqfile -t treefile -r rootfile -l splitfile`. There are no optional or extra flags; all four are needed for running the program.
`-s` : The input sequence file, in PHYLIP format.
`-t` : The input tree file, in Newick format. must match the taxa in the sequence file. Need not be rooted.
`-r` : The input root file for the tree, indicated by a list of taxa on one side of the root.
       Need not necessarily match the root indicated by the Newick tree, if it has one.
       Taxa are signified by integer numbers starting from 0, according to their order of apperance in the sequence file.
`-l` : The input split file indicating the FunDi split point by a list of taxa on one side of the split.
       Taxa are signified in the same way as `-r`.

- fundicalc
The command for fundicalc looks like `fundicalc -s seqfile -t treefile -l splitfile -m model -g gridsize`. All are required except for `-g`, which has a default value of `4`. 
`-s` : The input sequence file, in PHYLIP format.
`-t` : The name of the tree file that was used in rtsplit. While the file itself is not needed, the names of trees A, B, and AB are derived from the name of this file. 
`-l` : The input split file that was used in rtsplit. 
`-m` : The likelihood model to be used in IQ-TREE. E.g.: `LG+C60+G`
`-g` : The size of the grid to search for finding optimal split edge length `t`. DEFAULT: 4

### Output
`rtsplit` outputs three files: tree AB, tree A, and tree B. The name of these files is dependent on the name of the original tree file used in the program's input. NOTE: Do not rename these output files, as it may make fundicalc unusable.

`fundicalc` does not output any files as its answer is two numbers: `rho` and `t`. These values are indicated on the command line after the program has finished running. There are temporary files generated during runtime that begin with `tmp.`; if any are present after running the script they can be safely deleted.

### References

[fundi]

Susko, E. (2022). gfmix: Phylogenetic analyses using the site-and-branch-heterogeneous GFmix model Version 1.1: https://www.mathstat.dal.ca/~tsusko/doc/gfmix.pdf

Nguyen L.T., Schmidt H.A., von Haeseler A., Minh B.Q. (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274.
