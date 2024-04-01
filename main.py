import numpy as np
import pandas as pd
from os import walk
import re
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from os.path import join
import json
from concurrent.futures import ThreadPoolExecutor

def read_data_for_interaction(mypath: str, interaction_type: str) -> pd.DataFrame:
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
    data_list = []
    
    # Iterate over each filename
    for filename in filenames:
        # Extract DD from the filename
        numbers = re.findall(r'\d+',filename)
        HD = round(int(numbers[0]) / 42.0,2)
        DD = int(numbers[2]) / 1000.0
        variant_HD = int(numbers[1])
        variant_DD = int(numbers[3])
        position = int(numbers[4])
        file_path = join(mypath, filename)
        
        # Read data based on interaction type
        if interaction_type == 'Binding Energy':
            data = pd.read_table(file_path, delim_whitespace=True, header=0, usecols=[1], nrows=31).mean().round(2)
        elif interaction_type == 'Hydrogen Bonds':
            data = pd.read_table(file_path, delim_whitespace=True, header=70, usecols=list(range(27, 35)), nrows=31).mean().round(2)
        elif interaction_type == 'Hydrophobic Interactions':
            data = pd.read_table(file_path, delim_whitespace=True, header=70, usecols=list(range(35, 43)), nrows=31).mean().round(2)
        elif interaction_type == 'Ionic Interactions':
            data = pd.read_table(file_path, delim_whitespace=True, header=70, usecols=list(range(43, 51)), nrows=31).mean().round(2)
        else:
            raise ValueError("Invalid interaction type.")
            
        data_list.append((DD, HD, variant_HD, variant_DD, position, *data.values))
    
    # Convert data list to DataFrame
    if interaction_type == 'Binding Energy':
        df = pd.DataFrame(data_list, columns=['DD', 'HD', 'Variant_HD', 'Variant_DD', 'Position', 'Binding_Energy'])
    else:
        amino_acids = ["ALA", "ARG", "GLN", "GLU", "GLY", "HYP", "LEU", "PRO"]
        columns = ['DD', 'HD', 'Variant_HD', 'Variant_DD', 'Position'] + [f"{interaction_type}_{aa}" for aa in amino_acids]
        df = pd.DataFrame(data_list, columns=columns)

    return df

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



def process_directory(element, interaction_type):
    return read_data_for_interaction(element, interaction_type)

def hydroxylation_degree_analysis(HD: list, interaction_type: str) -> pd.DataFrame:
    """
    Extracting data for a given interaction type for selected collagen hydroxylation degrees
    """
    mypath = r"C:\Users\piotr\Downloads\Collagen+Chitosan\Collagen+Chitosan"
    subdir = "Bindenergy\\" if interaction_type == "Binding Energy" else "Analysis\\"
    
    filtered_directories = [join(mypath, f"Wyniki_{round(42 * h)}HYP", subdir) for h in HD]
    #print(filtered_directories)
    
    # Initialize a ThreadPoolExecutor with a maximum of 8 threads (adjust as needed)
    with ThreadPoolExecutor(max_workers=8) as executor:
        # Submit the processing of each directory to the executor
        futures = [executor.submit(process_directory, element, interaction_type) for element in filtered_directories]
        
        # Gather results as they become available
        data_list = [future.result() for future in futures]
    
    return pd.concat(data_list, ignore_index=True)


def main():
    HD = [0, 0.14, 0.29, 0.43, 0.57, 0.71, 0.86, 1]
    bind_energy = hydroxylation_degree_analysis(HD, "Binding Energy")
    hbond = hydroxylation_degree_analysis(HD, "Hydrogen Bonds")
    hydrophobic = hydroxylation_degree_analysis(HD, "Hydrophobic Interactions")
    ionic = hydroxylation_degree_analysis(HD, "Ionic Interactions")
    result = pd.merge(hbond, ionic, on=['HD','DD','Variant_HD', 'Variant_DD','Position'], how="left")
    result = pd.merge(result, hydrophobic, on=['HD','DD','Variant_HD', 'Variant_DD','Position'], how="left")
    result = pd.merge(result,  bind_energy, on=['HD','DD','Variant_HD', 'Variant_DD','Position'], how="outer")

    #print(result)
    result.to_csv('./Files/All_interactions.csv')

if __name__=="__main__":
    main()