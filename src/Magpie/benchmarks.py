import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

from importlib.resources import files


def pysr_benchmark():
    X = np.random.standard_normal(size=(130, 5))
    y = 2 * np.cos(X[:, 3]) + np.abs(X[:, 1]) ** 1.2
    print(X.shape, y.shape)
    return train_test_split(X, y, test_size=0.3, random_state=0)
    

def winequality_red():
    # use importlib to access a data file, even if the user
    # has installed using pip (ie doesn't have this csv in their current directory)
    fname = str(files(__package__) / 'data' / 'winequality-red.csv')
    d = pd.read_csv(fname, sep=';')
    X = d[d.columns[:-1]].values
    X = X / X.max(axis=0)
    y = d[d.columns[-1]].values
    y = y / y.max()
    return train_test_split(X, y, test_size=0.7, random_state=0)
