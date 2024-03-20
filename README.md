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
