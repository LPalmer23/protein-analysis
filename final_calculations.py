import pandas as pd  # Import pandas for data manipulation and analysis

# Load the Excel data file into a DataFrame
df = pd.read_excel('Schwanhäusser2011-SupplementaryTable3.xlsx')

# Define relevant columns in the DataFrame for easy access
protein_length_column = 'Protein length [amino acids]'
protein_copy_number_column = 'Protein copy number average [molecules/cell]'
mrna_copy_number_column = 'mRNA copy number average [molecules/cell]'
translation_rate_column = 'translation rate constant (ksp) average [molecules/(mRNA*h)]'

# Define the cell doubling time as a global variable
doubling_time_hours = 27  # Assumed cell doubling time in hours

# Function to calculate the mean protein size
def calculate_mean_protein_size(df):
    total_contribution = 0  # Initialize a variable to store the total contribution of protein sizes
    total_protein_copies = 0  # Initialize a variable to store the total number of protein copies

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Extract protein length and copy number from the current row
        protein_length = row[protein_length_column]
        protein_copy_number = row[protein_copy_number_column]
        
        # Ensure there are no missing values
        if pd.notna(protein_length) and pd.notna(protein_copy_number):
            # Calculate the contribution of the current protein
            contribution = protein_length * protein_copy_number
            total_contribution += contribution  # Add to total contribution
            total_protein_copies += protein_copy_number  # Add to total protein copies

    # Calculate the mean protein size if there are protein copies available
    if total_protein_copies > 0:
        mean_protein_size = total_contribution / total_protein_copies
        return mean_protein_size  # Return the calculated mean protein size
    else:
        return None  # Return None if there are no protein copies

# Step 2: Estimate Missing mRNA Copy Numbers

# Approach 1: Best Possible Translation Rate
def estimate_mrna_with_best_rate(df):
    best_translation_rate = 1 / 3  # Assume 1 protein synthesized every 3 seconds
    proteins_per_hour = 3600 * best_translation_rate  # Calculate proteins produced per hour
    
    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Extract protein and mRNA copy numbers
        protein_copy_number = row[protein_copy_number_column]
        mrna_copy_number = row[mrna_copy_number_column]
        
        # If mRNA copy number is missing, estimate it
        if pd.isna(mrna_copy_number):
            # Calculate the maximum number of mRNA needed to produce proteins
            max_mrna_needed = protein_copy_number / (proteins_per_hour * doubling_time_hours)
            # Update the DataFrame with the estimated mRNA copy number
            df.at[index, mrna_copy_number_column] = max_mrna_needed
    
    return df  # Return the updated DataFrame

# Approach 2: Median Translation Rate
def estimate_mrna_with_median_rate(df):
    # Filter out unrealistic translation rates and compute median
    valid_rates = df[translation_rate_column].dropna()  # Remove NaN values
    valid_rates = valid_rates[valid_rates <= 3600 / 2]  # Filter rates above the theoretical maximum
    median_translation_rate = valid_rates.median()  # Calculate the median translation rate

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Extract protein and mRNA copy numbers
        protein_copy_number = row[protein_copy_number_column]
        mrna_copy_number = row[mrna_copy_number_column]
        
        # If mRNA copy number is missing, estimate it
        if pd.isna(mrna_copy_number):
            # Calculate mRNA needed using the median translation rate
            median_mrna_needed = protein_copy_number / (median_translation_rate * doubling_time_hours)
            # Update the DataFrame with the estimated mRNA copy number
            df.at[index, mrna_copy_number_column] = median_mrna_needed
    
    return df  # Return the updated DataFrame

# Step 3: Normalize by Translation Rates
def normalize_mrna_by_translation_rate(df):
    reference_translation_rate = df[translation_rate_column].median()  # Determine reference translation rate

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Extract the translation rate
        translation_rate = row[translation_rate_column]
        # Normalize mRNA copy numbers if translation rate is valid
        if pd.notna(translation_rate) and translation_rate > 0:
            # Calculate normalization factor
            normalization_factor = reference_translation_rate / translation_rate
            # Adjust mRNA copy number by normalization factor
            df.at[index, mrna_copy_number_column] *= normalization_factor
    
    return df  # Return the normalized DataFrame

# Function to calculate the median mRNA coding sequence length
def calculate_median_mrna_length(df):
    mrna_lengths = []  # List to store mRNA lengths

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        # Extract mRNA length (using protein length as proxy) and mRNA copy number
        mrna_length = row[protein_length_column]
        mrna_copy_number = row[mrna_copy_number_column]

        # Append mRNA lengths based on copy number
        if pd.notna(mrna_length) and pd.notna(mrna_copy_number):
            for _ in range(int(mrna_copy_number)):
                mrna_lengths.append(mrna_length)

    mrna_lengths.sort()  # Sort mRNA lengths to prepare for median calculation
    n = len(mrna_lengths)  # Determine the number of mRNA lengths

    # Calculate and return the median mRNA length
    if n > 0:
        if n % 2 == 1:
            median_mrna_length = mrna_lengths[n // 2]
        else:
            mid_index = n // 2
            median_mrna_length = (mrna_lengths[mid_index - 1] + mrna_lengths[mid_index]) / 2
        return median_mrna_length  # Return the calculated median
    else:
        return None  # Return None if no mRNA lengths are available

# Main Execution
if __name__ == "__main__":
    # Calculate Mean Protein Size
    mean_protein_size = calculate_mean_protein_size(df)
    print("Mean Protein Size:", mean_protein_size)  # Output the mean protein size

    # Estimate mRNA with Best Translation Rate
    df_best_rate = estimate_mrna_with_best_rate(df.copy())  # Use the best rate approach
    median_mrna_length_best_rate = calculate_median_mrna_length(df_best_rate)  # Calculate median mRNA length
    print("Median mRNA Coding Sequence Length (Best Rate):", median_mrna_length_best_rate)

    # Estimate mRNA with Median Translation Rate
    df_median_rate = estimate_mrna_with_median_rate(df.copy())  # Use the median rate approach
    median_mrna_length_median_rate = calculate_median_mrna_length(df_median_rate)  # Calculate median mRNA length
    print("Median mRNA Coding Sequence Length (Median Rate):", median_mrna_length_median_rate)

    # Normalize by Translation Rates
    df_normalized = normalize_mrna_by_translation_rate(df.copy())  # Normalize mRNA data
    normalized_median_mrna_length = calculate_median_mrna_length(df_normalized)  # Calculate normalized median
    print("Normalized Median mRNA Coding Sequence Length:", normalized_median_mrna_length)