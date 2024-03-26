import numpy as np
import pandas as pd
from os import walk
import re
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from os.path import join
import json

def read_data_for_interation(mypath: str,interaction_type: str) -> dict:
    """
    Data for given interaction extraction, i.e. second column from *bindenergy_Mg.tab files for binding
    energies, for other interactions we take 8 columns that represent contribution of each amino acid
    building blocks of collagen variant. We look specifically into HD which is a #HYP/(#HYP + #PRO) amino
    acids. Note: there are ~1500 files depending on the case. The values from simulations
    take following range of HD values: [0.00, 0.14, 0.29, 0.43, 0.57, 0.71, 0.86, 1.00]. 
    There are three cases where the is one variant of collagen HD: 0(no hydroxylation), 0.43(native) 
    and 1(full hydroxylation). In other cases there are 5 variants of collagen for each HD. 
    In total there are around 40000 structures to be evaluated. There is binding energy
    file and analysis file. The latter contains information about hydrogen bonds, ionic and 
    hydrophobic interactions. Due to limitations on computations on supercomputers. The calculation of
    binding energy is much more demanding so this quantity was measured for much shorter. 
    Additionally the data can be gathered for DD for chitosan [0.125:1:0.125]
    These function extracts specific data for specific case and returns dictionary of all structures specified
    interaction.
    """
    filenames = next(walk(mypath), (None, None, []))[2]  # find all files from directory [] if no file
    pattern = '[1-9]\d{2,3}' # find values that repersent different deacetylation degrees 125:1000:125
    # Analyze deacetylation degree: match pattern below (add new variable to function) and goes over all 
    # directories Wyniki_* to extract specific DD 
    #pattern = rf'{pattern_variable}'
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
    """
    Drawing the data with best fitting line. Due to dense data set it is advised to plot max 4 lines
    """
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
        #equation = f'Y = {slope:.4f}X + {intercept:.2f}'
        # Add R2 value below the first point
        plt.text(x_values[6], 0.95 * y_means[6], r'$R^2 = {:.2f}$'.format(r_squared), fontsize=10, ha='center')

    plt.title(interaction_type + " vs deacetylation degree")
    plt.xlabel("Deacetylation degree [%]")
    plt.ylabel(interaction_type)
    plt.legend(loc='lower right')

    path = r"./Files/" + interaction_type.replace(' ','') + r"_vs_dd.png"
    plt.savefig(path)

def deacetylation_degree_analysis():
    """
    Extracting data for a given interaction type for selected chitosan deacetylation degrees
    """
    pass
def save_cleaned_data_to_file(data, filename, HD, interaction_type):
    """Save extracted data to JSON file for ML model"""
    AA = ['ALA', 'ARG', 'GLN', 'GLU', 'GLY', 'HYP', 'LEU', 'PRO']
    with open(filename, 'w') as json_file:
        all_data = {}
        for i, dictionary in enumerate(data):
            cleaned_data = {}
            for degree, values in dictionary.items():
                if interaction_type != 'Binding Energy':
                    # Use AA keys for all types of interactions except 'Binding Energy'
                    cleaned_values = [{'AA': {AA[j]: value.tolist() if isinstance(value, pd.Series) else value for j, (key, value) in enumerate(zip(AA, degree_values))}} for degree_values in values]
                else:
                    cleaned_values = []
                    for item in values:
                        # Convert any remaining Pandas Series objects to lists
                        cleaned_item = {key: value.tolist() if isinstance(value, pd.Series) else value for key, value in item.items()}
                        # Convert any nested Series to lists
                        cleaned_item = {key: [sub_item.tolist() if isinstance(sub_item, pd.Series) else sub_item for sub_item in value] if isinstance(value, list) else value for key, value in cleaned_item.items()}
                        cleaned_values.append(cleaned_item)
                    
                cleaned_data[degree] = cleaned_values
            
            # Add data for each HD value to the JSON structure
            all_data[str(HD[i])] = cleaned_data
        
        json.dump(all_data, json_file, indent=4)


def hydroxylation_degree_analysis(HD: list, interaction_type: str) -> list:
    """
    Extracting data for a given interaction type for selected colagen hydroxylation degrees
    """
    HD_cases = ["Wyniki_" + str(int(42 * x)) + "HYP" for x in HD] 
    mypath = r"C:\Users\piotr\Downloads\Collagen+Chitosan\Collagen+Chitosan"
    subdir = "Bindenergy\\" if interaction_type == "Binding Energy" else "Analysis\\"
    
    filtered_directories = [join(mypath, path, subdir) for path in HD_cases]
    data = []
    for element in filtered_directories:
        d = read_data_for_interation(element, interaction_type)
        data.append(d)
    return data

def amino_acids(a):
    """
    This function is assigning amino acid contribution to selected interaction. 
    Doesn't work for binding energy
    """
    AA = ['ALA' 'ARG' 'GLN' 'GLU' 'GLY' 'HYP' 'LEU' 'PRO']
   

def main():
    interaction_type = "Ionic Interactions"
    HD = [0, 0.14, 0.29, 0.43, 0.57, 0.71, 0.86, 1]
    f1 = hydroxylation_degree_analysis(HD, interaction_type)
    #create_figure(f1, HD, interaction_type)
    save_cleaned_data_to_file(f1,'./Files/Ionic_Interactions.json', HD, interaction_type)
if __name__=="__main__":
    main()
