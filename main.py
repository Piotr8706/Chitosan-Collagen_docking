import numpy as np
import pandas as pd


def read_data(file):
    data1 = pd.read_table(file, delim_whitespace=True)
    return data1

def main():
    file = r"C:\Users\piotr\Downloads\Collagen+Chitosan\Collagen+Chitosan\Wyniki_18HYP\Bindenergy\CollagentypeII18HYP1+Chitosan875_5_040_bindenergy_Mg.tab"


if __name__=="__main__":
    main()
