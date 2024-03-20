# MI-GO

## Introduction

MI-GO is a tool oriented to retrieve known GO-TERMS from gene2ids annotation files and then perform calculations like Information content (IC) or Shannon Diversity Index accross multiple datasets, which is useful for comparing GO-TERMs annotations across species and/or annotations and select genes and GO-TERMs based on if they are informative or not.

## Installation 

First, get MI-GO repository:

    git clone https://github.com/victorgcb1987/Mi-GO.git
    cd MI-GO

Then install the python requirements. A python venv is strongly advised.

    pip install -r requirements.txt
    python setup.py install

## Usage

A TAB-file with dataset's name and file path is needed in order to run the program. This file should be something like this:
    
    annot1    /path/to/annot1.txt
    annot2    /path/to/annot2.txt

Each dataset file should have the following format:

    AT1G12200.1|PACid:19651062	GO:0000325;GO:0050661;GO:0050660;GO:0050832;GO:0004499
    AT3G48680.1|PACid:19658684	GO:0010228;GO:0000325;GO:0009853;GO:0005515;GO:0009651;GO:0031966;GO:0046872;GO:0005747;GO:0005739;GO:0005829;GO:0009737
    AT1G14260.1|PACid:19654319	GO:0016021;GO:0000325;GO:0016740;GO:0046872;GO:0016567;GO:0016020

Also, a GO ontology file is needed, for example the go-basic.obo file found in https://geneontology.org/docs/download-ontology/

You can run MI-GO like this:

    python MI-GO.py -i {ListGOFiles.txt} -b {go-basic.obo.txt} -r -o {output_dir}

The following arguments are mandatory: -i/--input is the path to tab-delimited file with the list of filepaths to each dataset. -b/--obo is the filepath is to the GO ontology in .obo format file containing the description and status of each GO Term and -o/--output is directory path where all the results are going to be stored. -r/--remove obsolete is an optional argument which will remove from final results each GO Term obosolete by obo file's standard. 

## Output

After running the program, results should be stored in the {outp_dir} defined in the command line. If {output_dir} doesn't exist, the directory is created. Inside two types of files can be found:

* Diversity_IC_table.tsv: a tab-delimited file containing results for all species. It should look like this:

| NS | GO_ID      | dcnt | Depth | SP1_IC           | SP2_IC           | SP3_IC           | Diversity_IC     | Specifity_IC        | SP1_GeneCount | SP2_GeneCount | SP3_GeneCount | Diversity_GeneCount | Specifity_GeneCount | Name                         |
|----|------------|------|-------|------------------|------------------|------------------|------------------|---------------------|---------------|---------------|---------------|---------------------|---------------------|-----------------------------|
| BP | GO:0000003 | 31   | 1     | 5.613 | 6.184  | 5.774 | 1.097 | 0.001 | 21            | 19            | 31            | 1.074    | 0.034  | reproduction                |
| BP | GO:0002376 | 566  | 1     | 3.902 | 3.627 | 3.566   | 1.097 | 0.001 | 7             | 5             | 4             | 1.071    | 0.038  | immune system process       |
| MF | GO:0003774 | 9    | 1     | 5.483 | 5.563 | 4.307 | 1.092 | 0.009 | 0             | 1             | 1             | 0.693   | 0.584   | cytoskeletal motor activity |

NS indicates is GO Term is a Biological process (BP), Molecular Function (MF) or Cellular Component (CC). ICs are the Information content for each GO Term and species. Next two columns are for Shannon Diversity Index and specifity for ICs pertaining to this GO Term and species. The same is for the following columns, but instead of using IC uses the number of Genes in the annotation with this GO term assigned. Last column, Name, is the description of the GO Term.

* {Species}_SUM_IC.tsv: It's a tab-delimited file for each species. For each gene, it contains the IC sum of all GOI terms that are annotated for that gene. It should look like this:

| gene  | Sum_IC             |
|-------|--------------------|
| Gene1 | 37.982  |
| Gene2 | 22.270 |
| Gene2 | 36.629  |
