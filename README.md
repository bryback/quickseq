quickseq
========
QuickSeq is a simple tool that facilitates the design of homologous recombination sequences containing specific restriction sites for Recombinase-mediated cassette exchange and gene deletion in microbes.

The user inputs a genome file and a list of genes/operons to be deleted (e.g. lacZ, araBAD). The program returns annotated genbank files containing DNA sequences that are homologous to defined upstream and downstream regions and contain specific restriction sites (see user guide for details).

The front-end is written in python and uses webapp2 and Jinja2 to pass values from HTML to the back-end. jQuery and AJAX are used to dynamically update the input form (e.g. fill the autocomplete gene list after parsing the genome file).

The back-end is written in python and represented as a class called Genescript. Genomes and DNA sequences are parsed and manipulated using the Biopython library. Output is generated via StringIO and zipfile.

An instance of the current version of QuickSeq is hosted at quickseq0.appspot.com 
