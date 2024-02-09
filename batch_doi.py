import subprocess
import os
import sys


def run_simulation():
    # Define the output folder
    output_folder = "output"

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define the parameter ranges
    algorithm_range = [0, 1]
    num_sensors_range = [50, 100, 200, 400]
    sensor_radius_doi_percentage_range = [0.9, 0.85, 0.9]
    doi_range = [0.01, 0.05, 0.1]

    # Define the base command
    baseCommand = ".\cmake-build-release\TOSN.exe --params"

    # Initialize a counter for total iterations
    it = 1
    tot = (len(algorithm_range) * len(num_sensors_range) * len(sensor_radius_doi_percentage_range) * len(doi_range))

    for num_sensors in num_sensors_range:
        for algorithm in algorithm_range:
            for sensor_radius_doi_percentage in sensor_radius_doi_percentage_range:
                for doi in doi_range:
                    # Define the exp_name parameter based on the parameter values
                    expName = f"doi_s{num_sensors}_d1_r40_rd{sensor_radius_doi_percentage}_doi{doi}_b1.5_w1_a{algorithm}"

                    # Construct the full command with varying parameters
                    fullCommand = f"{baseCommand} -exp_name {expName} -iterations 10 -num_sensors {num_sensors} -num_depots 1 -sensor_radius 40 -energy_budget 1500000 -wireless_technology 1 -algorithm {algorithm} -seed {it} -scenario 0"

                    # Display the command
                    print(f"Executing {it}/{tot} = {fullCommand}")

                    # Uncomment the next line to execute the command
                    subprocess.run(fullCommand, shell=True)

                    # Increment the total iterations counter
                    it += 1


if __name__ == "__main__":
    run_simulation()
