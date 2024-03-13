from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from src.godag import (calculate_go_terms_IC_diversity,
                       calculate_go_terms_IC_specifity,
                       calculate_go_terms_geneCount_diversity,
                       calculate_go_terms_geneCount_specifity,
                       calculate_sum_IC,
                       convert_counts_stats_to_data_frame,
                       count_genes_by_go,
                       get_annotations,
                       get_subdag,
                       get_subdag_statistics,
                       get_terms_counts,
                       group_genes_by_GO,
                       read_godag,
                       write_sum_IC_tables)

def parse_arguments():
    desc = "MI-GO: tool for calculate IC content and Diversity for GO TERMS"
    parser = ArgumentParser(description=desc)

    help_input_fpath = "(Required) Path to a TAB-file with species and gene2go files paths"
    parser.add_argument("--input", 
                        "-i", type=str,
                        help=help_input_fpath,
                        required=True)
    help_obo_fpath = "(Required) Path to a go-basic in .obo format"
    parser.add_argument("--obo", 
                        "-b", type=str,
                        help=help_obo_fpath,
                        required=True)
    help_output_fpath = "(Required) Output Dir"
    parser.add_argument("--output", 
                        "-o", type=str,
                        help=help_output_fpath,
                        required=True) 
    return parser


def get_arguments():
    options = parse_arguments().parse_args()
    datasets = {}
    with open(options.input) as input_fhand:
        for line in input_fhand:
            if line:
                line = line.rstrip().split()
                datasets[line[0]] = Path(line[1])
        out_path = Path(options.output)
        if not out_path.exists():
            out_path.mkdir(parents=True)
    obo = Path(options.obo)
    return {"datasets": datasets,
            "obo": obo,
            "output": out_path} 
     
def main():
    arguments = get_arguments()
    godag = read_godag(arguments["obo"])
    gene_groups = group_genes_by_GO(arguments["datasets"])
    gene_counts = count_genes_by_go(gene_groups)
    annotations = get_annotations(godag, arguments["datasets"])
    counts = get_terms_counts(godag, datasets=annotations)
    gosubdag = get_subdag([goid for goid in godag], godag)
    sorted_nts = get_subdag_statistics(gosubdag, [go for go in godag], sort_by_depth=True)
    dataframe = convert_counts_stats_to_data_frame(sorted_nts, counts, gene_counts)
    calculate_go_terms_IC_diversity(dataframe)
    calculate_go_terms_IC_specifity(dataframe)
    calculate_go_terms_geneCount_diversity(dataframe)
    calculate_go_terms_geneCount_specifity(dataframe)
    sum_data = calculate_sum_IC(dataframe, annotations)
    write_sum_IC_tables(sum_data, arguments["output"])
    data_frame_fpath = arguments["output"] / "Diversity_IC_table.tsv"
    dataframe.to_csv(data_frame_fpath, sep="\t", index=False)

if __name__ == "__main__":
    main()