import subprocess
import os
import sys


def run_simulation(wireless_technology_value):
    # Define the output folder
    output_folder = "output"

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    baseline = False

    # Define the parameter ranges
    algorithm_range = [0, 1, 2, 3]
    num_sensors_range = [50, 100]
    num_depots_range = [1, 3, 5]
    sensor_radius_range = [20, 40, 60, 80]
    energy_budget_range = [1500000, 2000000, 2500000]
    wireless_technology_range = [wireless_technology_value]

    if baseline:
        sensor_radius_range = [0]
        algorithm_range = [0, 2]

    # Define the base command
    baseCommand = "./cmake-build-release/TOSN --params"

    # Initialize a counter for total iterations
    it = 1
    tot = (len(wireless_technology_range) * len(algorithm_range) * len(num_sensors_range) * len(num_depots_range) * len(sensor_radius_range) * len(energy_budget_range)) // 2

    for num_sensors in num_sensors_range:
        for num_depots in num_depots_range:
            for sensor_radius in sensor_radius_range:
                for energy_budget in energy_budget_range:
                    for algorithm in algorithm_range:
                        for wireless_technology in wireless_technology_range:
                            # Check the constraints on numDepots and algorithm
                            if (num_depots == 1 and (algorithm == 0 or algorithm == 1)) or (num_depots > 1 and (algorithm == 2 or algorithm == 3)):
                                # Define the exp_name parameter based on the parameter values
                                expName = f"reg_s{num_sensors}_d{num_depots}_r{sensor_radius}_b{round(energy_budget / 1e+6, 2)}_w{wireless_technology}_a{algorithm}"

                                # Construct the full command with varying parameters
                                fullCommand = f"{baseCommand} -exp_name {expName} -iterations 10 -num_sensors {num_sensors} -num_depots {num_depots} -sensor_radius {sensor_radius} -energy_budget {energy_budget} -wireless_technology {wireless_technology} -algorithm {algorithm} -seed {it} -scenario 0"

                                # Display the command
                                print(f"Executing {it}/{tot} = {fullCommand}")

                                # Uncomment the next line to execute the command
                                subprocess.run(fullCommand, shell=True)

                                # Increment the total iterations counter
                                it += 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <wireless_technology_value>")
        sys.exit(1)

    wireless_technology_value = int(sys.argv[1])
    run_simulation(wireless_technology_value)
