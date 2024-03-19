import numpy as np
import pandas as pd
from os import walk
import re
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import time


def read_data_bindenergy(mypath: str,interaction_type: str) -> dict:
    """
    Binding energy extraction, i.e. second column from *bindenergy_Mg.tab files
    """
    filenames = next(walk(mypath), (None, None, []))[2]  # find all files from directory [] if no file
    pattern = '[1-9]\d{2,3}' # find values that repersent different deacetylation degrees 125:1000:125
    digits = np.unique([re.findall(pattern, filename)[0] for filename in filenames])
    result = {}
    # create dictionary of all binding energies from each deacetylation degree
    for digit in digits:
        files_with_digit = [filename for filename in filenames if digit in filename]
        bindenergies = []
        for filename in files_with_digit:
            match interaction_type:
                case 'Binding Energy':
                    data = pd.read_table(mypath + filename, delim_whitespace=True, header=0, usecols=[1], nrows=31).mean()
                case 'Hydrogen Bonds':
                    data = pd.read_table(mypath + filename, delim_whitespace=True, header=70, usecols=list(range(27, 35)), nrows=31).mean()
                case 'Hydrophobic Interactions':
                    data = pd.read_table(mypath + filename, delim_whitespace=True, header=70, usecols=list(range(35, 43)), nrows=31).mean()
                case 'Ionic Interactions':
                    data = pd.read_table(mypath + filename, delim_whitespace=True, header=70, usecols=list(range(43, 51)), nrows=31).mean()
            bindenergies.append(data)
        result[digit] = bindenergies

    return result

def create_figure(data: dict,interaction_type: str) -> None:
    # Convert dictionary keys to integers
    x_values = np.array([int(key)/10.0 for key in data.keys()])
    if interaction_type == 'Binding Energy':
        num = -1
    else:
        num = 8
    # Calculate average and standard deviation for each set of values
    y_means = np.array([num*np.mean(values) for values in data.values()])
    y_sems = np.array([np.std(values) / np.sqrt(len(values)) for values in data.values()])

    # Fit a linear regression line
    slope, intercept = np.polyfit(x_values, y_means, 1)

    # Calculate predicted y values
    y_pred = slope * x_values + intercept

    # Calculate R2 value
    r_squared = r2_score(y_means, y_pred)

    # Plot data points with error bars and the fitting line
    plt.figure()
    plt.errorbar(x_values, y_means, yerr=y_sems, fmt='o', label='Data Points with Error Bars')
    plt.plot(x_values, y_pred, color='red', label='Fitting Line')
    plt.title(interaction_type + " vs deacetylation degree")
    plt.xlabel("Deacetylation degree [%]")
    plt.ylabel(interaction_type)
    plt.legend(loc='lower right')

    # Add equation and R2 value to the plot
    equation = f'Y = {slope:.2f}X + {intercept:.2f}'
    r_squared_text = f'R2 = {r_squared:.2f}'
    plt.text(0.5, 0.95, equation, ha='center', va='center', transform=plt.gca().transAxes, fontsize=12)
    plt.text(0.5, 0.90, r_squared_text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=12)
    path = r"./Files/" + interaction_type + r"_vs_dd.png"
    plt.savefig(path)
    #plt.show()

def deacetylation_degree_analysis():
    pass

def hydroxylation_degree_analysis():
    pass

def amino_acids(a):
    AA = ['ALA' 'ARG' 'GLN' 'GLU' 'GLY' 'HYP' 'LEU' 'PRO']
   

def main():

    mypath_be = r"C:\Users\piotr\Downloads\Collagen+Chitosan\Collagen+Chitosan\Wyniki_36HYP_5\Bindenergy\\"
    mypath_inter = r"C:\Users\piotr\Downloads\Collagen+Chitosan\Collagen+Chitosan\Wyniki_42HYP\Analysis\\"
    start = time.time()
    result = read_data_bindenergy(mypath_inter,'Hydrogen Bonds')
    create_figure(result,'Hydrogen Bonds')
    #print(result)
    end = time.time()
    print((end-start))

if __name__=="__main__":
    main()
