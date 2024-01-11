# Define the parameter ranges
$algorithmRange = 0, 1, 2, 3
$numSensorsRange = 50, 100, 200, 400
$numDepotsRange = 1, 3, 5
$sensorRadiusRange = 20, 60, 80
$energyBudgetRange = 1500000, 2000000, 2500000
$wirelessTechnologyRange = 3

# Define the base command
$baseCommand = ".\cmake-build-release\TOSN.exe --params"

# Initialize a counter for total iterations
$it = 1
$tot = ($wirelessTechnologyRange.Count * $algorithmRange.Count * $numSensorsRange.Count * $numDepotsRange.Count * $sensorRadiusRange.Count * $energyBudgetRange.Count) / 2

# Nested loops to iterate over parameter combinations
foreach ($numSensors in $numSensorsRange) {
    foreach ($numDepots in $numDepotsRange) {
        foreach ($sensorRadius in $sensorRadiusRange) {
            foreach ($energyBudget in $energyBudgetRange) {
                foreach ($algorithm in $algorithmRange) {
					foreach ($wirelessTechnology in $wirelessTechnologyRange) {

						# Check the constraints on numDepots and algorithm
						if (($numDepots -eq 1 -and ($algorithm -eq 0 -or $algorithm -eq 1)) -or
							($numDepots -gt 1 -and ($algorithm -eq 2 -or $algorithm -eq 3))) {

							# Define the exp_name parameter based on the parameter values
							$expName = "reg_s${numSensors}_d${numDepots}_r${sensorRadius}_b$([math]::Round($energyBudget/1e+6,2))_w${wirelessTechnology}_a${algorithm}"

							# Construct the full command with varying parameters
							$fullCommand = "$baseCommand -exp_name $expName -num_sensors $numSensors -num_depots $numDepots -sensor_radius $sensorRadius -energy_budget $energyBudget -wireless_technology $wirelessTechnology -algorithm $algorithm -scenario 0"

							# Display the command
							Write-Host "Executing $it/$tot = $fullCommand"

							# Uncomment the next line to execute the command
							Invoke-Expression $fullCommand

							# Increment the total iterations counter
							$it++
						}
					}
                }
            }
        }
    }
}
