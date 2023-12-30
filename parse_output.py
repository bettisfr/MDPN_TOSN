import os
import pandas as pd


def merge_csv_files(output_file_path):
    folder_path='output'

    # Get a list of all CSV files in the directory starting with "reg_", "doi", or "dtr_"
    csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv') and file.startswith(('reg_', 'doi', 'dtr_'))]

    # Read the first CSV file to get the header
    first_file_path = os.path.join(folder_path, csv_files[0])
    header = pd.read_csv(first_file_path, nrows=0).columns.tolist()

    # Initialize an empty DataFrame to store the merged data
    merged_data = pd.DataFrame(columns=header)

    # Iterate through each CSV file and append its data to the merged_data DataFrame
    for csv_file in csv_files:
        file_path = os.path.join(folder_path, csv_file)

        # Read the CSV file into a DataFrame
        data = pd.read_csv(file_path)

        # Append the data to the merged_data DataFrame
        merged_data = pd.concat([merged_data, data], ignore_index=True)

    # Reset the index
    merged_data.reset_index(drop=True, inplace=True)

    # Create a new CSV file with the merged data
    merged_data.to_csv(output_file_path, index=False)

    print(f'Merged data saved to {output_file_path}')


def filter_and_save_subtable(input_file_path, output_file_path_prefix, scenario, energy_budget, sensor_radius, num_depots):
    # Read the original CSV file
    df = pd.read_csv(input_file_path)

    # Filter the data based on conditions
    filtered_df = df[
        (df['scenario'] == scenario) &
        (df['energy_budget'] == energy_budget) &
        (df['sensor_radius'] == sensor_radius) &
        (df['num_depots'] == num_depots)
        ]

    # Iterate over unique algorithm values and save a separate file for each
    for algorithm_value in filtered_df['algorithm'].unique():
        # Extract the suffix from the algorithm column
        suffix = f'_{algorithm_value}'

        # Filter data for the specific algorithm
        algorithm_df = filtered_df[filtered_df['algorithm'] == algorithm_value]

        # Sort the DataFrame by 'num_sensors' in ascending order
        algorithm_df_sorted = algorithm_df.sort_values(by='num_sensors', ascending=True)

        # Construct the output file path with suffix
        output_file_path = f'{output_file_path_prefix}{suffix}.csv'

        # Save the sorted DataFrame to CSV
        algorithm_df_sorted.to_csv(output_file_path, index=False)

        print(f'Subtable data for algorithm {algorithm_value} saved to {output_file_path}')


if __name__ == "__main__":
    file_path = 'merged_output.csv'

    # Define your specific conditions
    scenario_value = 0
    energy_budget_value = 2.5
    sensor_radius_value = 50
    num_depots_value = 1

    # Create subtables for each algorithm
    filter_and_save_subtable(file_path, 'subtable_output', scenario_value, energy_budget_value, sensor_radius_value, num_depots_value)
