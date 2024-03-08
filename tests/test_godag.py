import unittest

from pathlib import Path

from src.godag import (read_godag,
                       select_goids_from_godag,
                       get_annotations)


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
        print(annotations)
        