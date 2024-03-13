import pandas as pd

from math import log as ln
from math import     log2 as log

from goatools.anno.idtogos_reader import IdToGosReader
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.semantic import get_info_content, TermCounts


def convert_counts_stats_to_data_frame(dag_stats, counts, gene_counts):
    data_dict = {"NS": [], "GO_ID": [], "dcnt": [], "Depth": []}
    for species in counts:
        column_name = "{}_IC".format(species)
        data_dict[column_name] = []
    for species in counts:
        column_name = "{}_GeneCount".format(species)
        data_dict[column_name] = []
    data_dict["Name"] = []
    for values in dag_stats:
            data_dict["NS"].append(values.NS)
            data_dict["GO_ID"].append(values.GO)
            data_dict["dcnt"].append(values.dcnt)
            data_dict["Depth"].append(values.depth)
            for species, count in counts.items():
                column_name = "{}_IC".format(species)
                data_dict[column_name].append(get_info_content(values.GO, count))
            for species, count in counts.items():
                column_name = "{}_GeneCount".format(species)
                data_dict[column_name].append(gene_counts[species].get(values.GO, 0))
            data_dict["Name"].append(values.GO_name)
    dataframe = pd.DataFrame.from_dict(data_dict)
    return dataframe


def group_genes_by_GO(data_sets):
    go_groups = {species:{} for species in data_sets}
    for species, dataset in data_sets.items():
        with open(dataset) as fhand:
            for line in fhand:
                line = line.rstrip().split()
                gene = line[0]
                go_terms = line[1].split(";")
                for go_term in go_terms:
                    if go_term not in go_groups[species]:
                        go_groups[species][go_term] = [gene]
                    else:
                        go_groups[species][go_term].append(gene)
    return go_groups


def count_genes_by_go(grouped_genes):
    go_counts = {species: {} for species in grouped_genes}
    for species, go_terms in grouped_genes.items():
        for goterm, genes in go_terms.items():
            go_counts[species][goterm] = len(genes)
    return go_counts
                    

def read_godag(obo_fpath):
    godag = get_godag(obo_fpath)
    return godag


def select_goids_from_godag(godag, criteria=[], fields=[]):
    go_ids = []
    for go in godag.values():
        for criterion in criteria:
            for field in fields:
                if field == "name":
                    check = go.name
                if field == "namespace":
                    check = go.namespace
            if criterion in check:
                go_ids.append(go.item_id)
    return set(go_ids)


def get_godag_subset(ids, godag):
    return GoSubDag(ids, godag)


def get_annotations(godag, datasets={}):
    return {species : IdToGosReader(fpath, godag=godag, taxid=False) for species, fpath in datasets.items()}


def get_terms_counts(godag, datasets={}):
    return {species: TermCounts(godag, annot.get_id2gos_nss()) for species, annot in datasets.items()}


def get_subdag(go_list, godag):
    return GoSubDag(go_list, godag)


def get_subdag_statistics(subdag, go_terms, sort_by_depth=True):
    stats = [subdag.go2nt[go] for go in go_terms]
    if sort_by_depth:
        return sorted(stats, key=lambda nt: nt.depth)
    else:
        return stats


def calculate_go_terms_IC_diversity(dataframe):
    diversity = {"Diversity_IC": []}
    col_names = [colname for colname in dataframe.columns if "_IC" in colname and "Specifity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        values = [(float(raw_value)/N) * ln(float(raw_value)/N) if raw_value > 0 else 0 for raw_value in raw_values]
        diversity_value =  -sum(value for value in values if value != 0)
        diversity["Diversity_IC"].append(diversity_value)
    dataframe.insert(dataframe.columns.get_loc(col_names[-1])+1, 
                     "Diversity_IC", diversity["Diversity_IC"])
    return dataframe

def calculate_go_terms_IC_specifity(dataframe):
    specifity = {"Specifity_IC": []}
    col_names = [colname for colname in dataframe.columns if "_IC" in colname and "Diversity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        pijs = [float(raw_value/N) if N > 0 else 0 for raw_value in raw_values ]
        pi = float((1/len(raw_values))) * sum(pijs)
        values = [(pij/pi) * log(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
        si = (1/len(raw_values)) * sum(values)
        specifity["Specifity_IC"].append(si)
    dataframe.insert(dataframe.columns.get_loc("Diversity_IC")+1, 
                     "Specifity_IC", specifity["Specifity_IC"])
    return dataframe


def calculate_go_terms_geneCount_diversity(dataframe):
    diversity = {"_GeneCount": []}
    col_names = [colname for colname in dataframe.columns if "_GeneCount" in colname and "Specifity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        values = [(float(raw_value)/N) * ln(float(raw_value)/N) if raw_value > 0 else 0 for raw_value in raw_values]
        diversity_value =  -sum(value for value in values if value != 0)
        diversity["_GeneCount"].append(diversity_value)
    dataframe.insert(dataframe.columns.get_loc(col_names[-1])+1, 
                     "Diversity_GeneCount", diversity["_GeneCount"])
    return dataframe


def calculate_go_terms_geneCount_specifity(dataframe):
    specifity = {"Specifity_GeneCount": []}
    col_names = [colname for colname in dataframe.columns if "GeneCount" in colname and "Diversity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
        pi = float((1/len(raw_values))) * sum(pijs)
        values = [(pij/pi) * log(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
        si = (1/len(raw_values)) * sum(values)
        specifity["Specifity_GeneCount"].append(si)
    dataframe.insert(dataframe.columns.get_loc("Diversity_GeneCount")+1, 
                     "Specifity_GeneCount", specifity["Specifity_GeneCount"])
    return dataframe


def calculate_sum_IC(dataframe, annotations):
    sum_data = {species: {"gene": [], "Sum_IC":[]} for species in annotations}
    for species, annotation in annotations.items():
        for gene, goset in annotation.id2gos.items():
            species_ic_label = "{}_IC".format(species)
            sub_df = dataframe.get(["GO_ID", species_ic_label])
            sub_df = dict(zip(sub_df["GO_ID"], sub_df[species_ic_label]))
            if len(goset) == 0:
                sum_data[species]["gene"].append(gene)
                sum_data[species]["Sum_IC"].append(0)
            else:
                ic_sum = sum([sub_df[go_term] for go_term in goset])
                sum_data[species]["gene"].append(gene)
                sum_data[species]["Sum_IC"].append(ic_sum)
    return sum_data


def write_sum_IC_tables(sums_IC, out_fpath):
    for species, sum_IC in sums_IC.items():
        species_out_fpath = out_fpath / "{}_SUM_IC.tsv".format(species)
        df = pd.DataFrame.from_dict(sum_IC)
        df.to_csv(species_out_fpath, sep="\t", index=False)

