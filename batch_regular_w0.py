import subprocess
import os

# Define the output folder
output_folder = "output"

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# Define the parameter ranges
algorithmRange = [0, 1, 2, 3]
numSensorsRange = [50, 100, 200, 400]
numDepotsRange = [1, 3, 5]
sensorRadiusRange = [20, 60, 80]
energyBudgetRange = [1500000, 2000000, 2500000]
wirelessTechnologyRange = [0]

# Define the base command
baseCommand = ".\cmake-build-release\TOSN.exe --params"

# Initialize a counter for total iterations
it = 1
tot = (len(wirelessTechnologyRange) * len(algorithmRange) * len(numSensorsRange) *
       len(numDepotsRange) * len(sensorRadiusRange) * len(energyBudgetRange)) // 2

for numSensors in numSensorsRange:
    for numDepots in numDepotsRange:
        for sensorRadius in sensorRadiusRange:
            for energyBudget in energyBudgetRange:
                for algorithm in algorithmRange:
                    for wirelessTechnology in wirelessTechnologyRange:
                        # Check the constraints on numDepots and algorithm
                        if (numDepots == 1 and (algorithm == 0 or algorithm == 1)) or \
                                (numDepots > 1 and (algorithm == 2 or algorithm == 3)):

                            # Define the exp_name parameter based on the parameter values
                            expName = f"reg_s{numSensors}_d{numDepots}_r{sensorRadius}_b{round(energyBudget/1e+6, 2)}_w{wirelessTechnology}_a{algorithm}"

                            # Construct the full command with varying parameters
                            fullCommand = f"{baseCommand} -exp_name {expName} -iterations 10 -num_sensors {numSensors} -num_depots {numDepots} -sensor_radius {sensorRadius} -energy_budget {energyBudget} -wireless_technology {wirelessTechnology} -algorithm {algorithm} -seed {it} -scenario 0"

                            # Check if the CSV file for the experiment already exists
                            csv_file_path = os.path.join(output_folder, f"{expName}.csv")
                            if os.path.exists(csv_file_path):
                                # print(f"Skipping {it}/{tot} = {fullCommand}. Experiment already exists.")
                                a = 1
                            else:
                                # Display the command
                                print(f"Executing {it}/{tot} = {fullCommand}")

                                # Uncomment the next line to execute the command
                                subprocess.run(fullCommand, shell=True)

                            # Increment the total iterations counter
                            it += 1