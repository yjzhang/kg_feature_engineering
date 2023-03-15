# TODO: get info on available graphs
import os
import sys

import pandas as pd

PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(PATH, 'processed_graphs')

def get_available_graphs():
    files = os.listdir(DATA_PATH)
    return files


def load_graph(filename):
    files = os.listdir(DATA_PATH)
    if filename not in files:
        # if filename not in files, try opening it as a path
        if not os.path.exists(filename):
            raise FileNotFoundError()
        f = open(filename)
    else:
        f = open(os.path.join(PATH, filename))
    df = pd.read_csv(f)
    f.close()
    return df


def df_to_networkx(df):
    pass
