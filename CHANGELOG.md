# v2.0.0-beta
1. arguments of src/phylDiag.py and src/phylDiagViewer.py have changed. New names of arguments correspond to names in *High precision detection of conserved segments from synteny blocks*
2. arguments are now parsed with argparse instead of utils.myTools.checkArgs, of LibsDyogen
3. previous versions, with old arguments, can still be accessed via deprecated phylDiagCA.py and phylDiagViewerCA.py
4. new convention for extensions of files:
	* \*.genome.bz2 is a genome compressed with bz2 (instead of genes.\*.list.bz2)
	* \*.families.bz2 is set of families (instead of ancGenes.\*.list.bz2)
	* \*.sbs are synteny blocks, sometimes called \*.cs when they are conserved segments

# v2.0.0-alpha
