import unittest
import pandas as pd

from pathlib import Path

from src.godag import (convert_counts_stats_to_data_frame,
                       get_annotations,
                       get_subdag,
                       get_subdag_statistics,
                       get_terms_counts,
                       read_godag,
                       select_goids_from_godag)


class TestGodag(unittest.TestCase):


    def setUp(self):
        self.test_path = Path(__file__).parent.absolute() / "test_data"
    

    def test_select_goids_from_godag(self):
        obo_fpath = self.test_path / "go-basic.obo"
        godag = read_godag(obo_fpath)
        print(len([value.item_id for value in godag.values()]))
        fields=["namespace"]
        criteria = ["cellular_component", "molecular_function", "biological_process"]
        selected_ids = select_goids_from_godag(godag, criteria=criteria, fields=fields)
        print(len(selected_ids))
        criteria.pop(0)
        selected_ids = select_goids_from_godag(godag, criteria=criteria, fields=fields)
        print(len(selected_ids))
        criteria.pop(0)
        selected_ids = select_goids_from_godag(godag, criteria=criteria, fields=fields)
        print(len(selected_ids))
        criteria.pop()
        selected_ids = select_goids_from_godag(godag, criteria=criteria, fields=fields)
        print(len(selected_ids))


    def test_get_information_content(self):
        obo_fpath = self.test_path / "go-basic.obo"
        godag = read_godag(obo_fpath)
        datasets = {"human": self.test_path / "human.gene2go",
                    "mouse": self.test_path / "mouse.gene2go",
                    "fly": self.test_path / "fly.gene2go"}
        annotations = get_annotations(godag, datasets=datasets)
        counts = get_terms_counts(godag, datasets=annotations)
        gosubdag = get_subdag([goid for goid in godag], godag)
        sorted_nts = get_subdag_statistics(gosubdag, [go for go in godag], sort_by_depth=True)
        dataframe = convert_counts_stats_to_data_frame(sorted_nts, counts)
        print(dataframe)

        