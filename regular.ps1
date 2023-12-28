# Define the parameter ranges
$numSensorsRange = 100, 200, 300, 400, 500
$numDepotsRange = 1..5
$sensorRadiusRange = 50, 60, 70, 80, 90, 100
$energyBudgetRange = 2000000, 4000000, 6000000, 8000000, 10000000
$algorithmRange = 0..3

# Define the base command
$baseCommand = ".\cmake-build-release\TOSN.exe --params"

# Nested loops to iterate over parameter combinations
foreach ($numSensors in $numSensorsRange) {
    foreach ($numDepots in $numDepotsRange) {
        foreach ($sensorRadius in $sensorRadiusRange) {
            foreach ($energyBudget in $energyBudgetRange) {
                foreach ($algorithm in $algorithmRange) {
                    # Define the exp_name parameter based on the parameter values
                    $expName = "reg_s${numSensors}_d${numDepots}_r${sensorRadius}_b$([math]::Round($energyBudget/1e+6,2))_a${algorithm}"

                    # Construct the full command with varying parameters
                    $fullCommand = "$baseCommand -exp_name $expName -num_sensors $numSensors -num_depots $numDepots -sensor_radius $sensorRadius -energy_budget $energyBudget -algorithm $algorithm -scenario 0"

                    # Display the command
                    Write-Host "Executing: $fullCommand"

                    # Uncomment the next line to execute the command
                    # Invoke-Expression $fullCommand
                }
            }
        }
    }
}
