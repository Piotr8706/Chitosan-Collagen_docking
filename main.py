import numpy as np
import pandas as pd
from os import walk
import re
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from os.path import join
import glob

def read_data_bindenergy(mypath: str,interaction_type: str) -> dict:
    """
    Binding energy extraction, i.e. second column from *bindenergy_Mg.tab files.
    Note: there are ~1500 files depending on the case. There are three cases where the is one variant of 
    collagen HD: 0(no hydroxylation), 0.43(native) and 1(full hydroxylation). In other cases there are 5 variants
    of collagen for each HD. In total there are around 40000 structures to be evaluated. There is binding energy
    file and analysis file. The latter contains information about hydrogen bonds, ionic and 
    hydrophobic interactions. Due to limitations on computations on supercomputers. The calculation of
    binding energy is much more demanding so this quantity was measured for much shorter.

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

def create_figure(datasets: list, HD: list, interaction_type: str) -> None:
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'pink']  # Define colors for different datasets
    num_datasets = len(datasets)
    if num_datasets > len(colors):
        raise ValueError("Too many datasets. Please provide fewer datasets or add more colors.")

    plt.figure()

    for i, data in enumerate(datasets):
        # Convert dictionary keys to integers
        x_values = np.array([int(key)/10.0 for key in data.keys()])

        if interaction_type == 'Binding Energy':
            num = -1
        else:
            num = 8

        # Calculate average and standard deviation for each set of values
        y_means = np.array([num*np.mean(values) for values in data.values()])
        y_sems = np.array([np.std(values) / np.sqrt(len(values)) for values in data.values()])

        # Fit linear regression line
        slope, intercept = np.polyfit(x_values, y_means, 1)

        # Calculate predicted y values
        y_pred = slope * x_values + intercept

        # Calculate R2 value
        r_squared = r2_score(y_means, y_pred)
        # Plot data points with error bars and the fitting line
        plt.errorbar(x_values, y_means, yerr=y_sems, fmt='o', label=f'HD = {HD[i]}', color=colors[i])
        plt.plot(x_values, y_pred, color=colors[i], linestyle='-')

        # Add equation and R2 value to the plot
        equation = f'Y = {slope:.2f}X + {intercept:.2f}'
        r_squared_text = f'R2 = {r_squared:.2f}'
        plt.text(0.5, 0.35 + 0.2*i, equation, ha='center', va='center', transform=plt.gca().transAxes, fontsize=10, color=colors[i])
        plt.text(0.5, 0.30 + 0.2*i, r_squared_text, ha='center', va='center', transform=plt.gca().transAxes, fontsize=10, color=colors[i])

    plt.title(interaction_type + " vs deacetylation degree")
    plt.xlabel("Deacetylation degree [%]")
    plt.ylabel(interaction_type)
    plt.legend(loc='lower right')

    path = r"./Files/" + interaction_type.replace(' ','') + r"_vs_dd.png"
    plt.savefig(path)

def deacetylation_degree_analysis():
    pass

def hydroxylation_degree_analysis(HD: list, interaction_type: str) -> list:
    """
    Extracting data for a given interaction type for selected hydroxylation degrees
    """
    HD_cases = ["Wyniki_" + str(int(42 * x)) + "HYP" for x in HD] 
    mypath = r"C:\Users\piotr\Downloads\Collagen+Chitosan\Collagen+Chitosan"
    subdir = "Bindenergy\\" if interaction_type == "Binding Energy" else "Analysis\\"
    
    filtered_directories = [join(mypath, path, subdir) for path in HD_cases]
    data = []
    for element in filtered_directories:
        d = read_data_bindenergy(element, interaction_type)
        data.append(d)
    return data

def amino_acids(a):
    AA = ['ALA' 'ARG' 'GLN' 'GLU' 'GLY' 'HYP' 'LEU' 'PRO']
   

def main():
    interaction_type = "Hydrogen Bonds"
    HD = [0, 0.43, 1]
    f1 = hydroxylation_degree_analysis(HD, interaction_type)
    create_figure(f1, HD, interaction_type)

if __name__=="__main__":
    main()
