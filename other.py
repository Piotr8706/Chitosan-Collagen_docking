import json
import pandas as pd

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