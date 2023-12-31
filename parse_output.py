import os
import pandas as pd


def merge_csv_files():
    folder_path = 'output'
    output_folder = 'plot'

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

    # Create the "plot" folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Create a new CSV file with the merged data inside the "plot" folder
    output_file_path = os.path.join(output_folder, 'merged_output.csv')
    merged_data.to_csv(output_file_path, index=False)

    print(f'Merged data saved to {output_file_path}')


def filter_plot_res(energy_budget, sensor_radius, num_depots):
    input_file_path = 'plot/merged_output.csv'
    prefix = 'res'
    df = pd.read_csv(input_file_path)

    # Filter the data based on conditions
    filtered_df = df[
        (df['scenario'] == 0) &
        (df['energy_budget'] == energy_budget) &
        (df['sensor_radius'] == sensor_radius) &
        (df['num_depots'] == num_depots)
        ]

    os.makedirs('plot', exist_ok=True)

    for algorithm_value in filtered_df['algorithm'].unique():
        suffix = f'_d{num_depots}_r{int(sensor_radius)}_b{energy_budget:.1f}_a{algorithm_value}'

        algorithm_df = filtered_df[filtered_df['algorithm'] == algorithm_value]
        algorithm_df_sorted = algorithm_df.sort_values(by='num_sensors', ascending=True)

        output_file_path = os.path.join('plot', f'{prefix}{suffix}.csv')

        algorithm_df_sorted.to_csv(output_file_path, index=False)

        print(f'Saved to {output_file_path}')


def filter_plot_doi(energy_budget, sensor_radius, sensor_radius_doi_percentage, doi, num_depots):
    input_file_path = 'plot/merged_output.csv'
    prefix = 'doi'
    df = pd.read_csv(input_file_path)

    # Filter the data based on conditions
    filtered_df = df[
        (df['scenario'] == 1) &
        (df['energy_budget'] == energy_budget) &
        (df['sensor_radius'] == sensor_radius) &
        (df['sensor_radius_doi_percentage'] == sensor_radius_doi_percentage) &
        (df['doi'] == doi) &
        (df['num_depots'] == num_depots)
        ]

    os.makedirs('plot', exist_ok=True)

    for algorithm_value in filtered_df['algorithm'].unique():
        suffix = f'_d{num_depots}_r{int(sensor_radius)}_rd{sensor_radius_doi_percentage:.3f}_doi{doi:.3f}_b{energy_budget:.1f}_a{algorithm_value}'

        algorithm_df = filtered_df[filtered_df['algorithm'] == algorithm_value]
        algorithm_df_sorted = algorithm_df.sort_values(by='num_sensors', ascending=True)

        output_file_path = os.path.join('plot', f'{prefix}{suffix}.csv')

        algorithm_df_sorted.to_csv(output_file_path, index=False)

        print(f'Saved to {output_file_path}')


if __name__ == "__main__":
    merge_csv_files()

    # RES
    energy_budgets = [2.5]
    sensor_radii = [50, 75, 100]
    num_depots_values = [1, 3, 5]

    for energy_budget in energy_budgets:
        for sensor_radius in sensor_radii:
            for num_depots in num_depots_values:
                filter_plot_res(energy_budget, sensor_radius, num_depots)


    # DOI
    # filter_plot_doi(5, 50, 0.8, 0.001, 1)
