import pandas as pd

from goatools.anno.idtogos_reader import IdToGosReader
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.semantic import get_info_content, TermCounts



def convert_counts_stats_to_data_frame(dag_stats, counts):
    data_dict = {"NS": [], "GO_ID": [], "dcnt": [], "Depth": []}
    for species in counts:
        data_dict[species] = []
    data_dict["Name"] = []
    for values in dag_stats:
            data_dict["NS"].append(values.NS)
            data_dict["GO_ID"].append(values.GO)
            data_dict["dcnt"].append(values.dcnt)
            data_dict["Depth"].append(values.depth)
            for species, count in counts.items():
                data_dict[species].append(get_info_content(values.GO, count))
            data_dict["Name"].append(values.GO_name)
    dataframe = pd.DataFrame.from_dict(data_dict)
    return dataframe


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

