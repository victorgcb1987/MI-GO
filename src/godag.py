from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.anno.idtogos_reader import IdToGosReader


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
    #return {species : get_objanno(fpath, anno_type="gene2go", godag=godag, taxid=False) for species, fpath in datasets.items()}
    return {species : IdToGosReader(fpath, godag=godag, taxid=False) for species, fpath in datasets.items()}


def calculate_information_content(godag, datasets={}, format="gene2go"):
    pass