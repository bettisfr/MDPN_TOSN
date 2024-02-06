#include "output.h"

// Function to calculate the average of a vector of doubles
double calculate_average(const vector<double> &values) {
    if (values.empty()) {
        return 0.0; // Handle the case where there are no values to average
    }

    double sum = 0.0;
    for (double value: values) {
        sum += value;
    }

    return sum / values.size();
}

// Function to calculate the standard deviation of a vector of doubles
double calculate_stddev(const vector<double> &values, double average) {
    if (values.size() <= 1) {
        return 0.0; // Handle the case where there is insufficient data for standard deviation
    }

    double sum_squared_diff = 0.0;
    for (double value: values) {
        double diff = value - average;
        sum_squared_diff += diff * diff;
    }

    return sqrt(sum_squared_diff / (values.size() - 1));
}

void save_output(const input &par, vector<solution> &results) {
    string filename = "output/" + par.exp_name + ".csv";
    ofstream file(filename, ofstream::out | ofstream::trunc);

    // Initialize averages and standard deviations
    double tours_number_avg = -1;
    double uncovered_sensors_avg = -1;
    double lost_data_avg = -1;
    double running_time_avg = -1;

    double tours_number_std = -1;
    double uncovered_sensors_std = -1;
    double lost_data_std = -1;
    double running_time_std = -1;

    // Calculate averages and standard deviations
    if (!results.empty()) {
        vector<double> tours_numbers;
        vector<double> uncovered_sensors;
        vector<double> lost_data;
        vector<double> running_times;

        for (const auto &out: results) {
            tours_numbers.push_back(out.tours_number);
            uncovered_sensors.push_back(out.uncovered_sensors * 100.0 / par.num_sensors); // Calculate as a percentage
            lost_data.push_back(out.lost_data * 100.0 / out.total_data); // Calculate as a percentage
            running_times.push_back(out.running_time);
        }

        tours_number_avg = calculate_average(tours_numbers);
        uncovered_sensors_avg = calculate_average(uncovered_sensors);
        lost_data_avg = calculate_average(lost_data);
        running_time_avg = calculate_average(running_times);

        tours_number_std = calculate_stddev(tours_numbers, tours_number_avg);
        uncovered_sensors_std = calculate_stddev(uncovered_sensors, uncovered_sensors_avg);
        lost_data_std = calculate_stddev(lost_data, lost_data_avg);
        running_time_std = calculate_stddev(running_times, running_time_avg);
    }

    if (file.is_open()) {
        file
                << "num_sensors,num_depots,sensor_radius,sensor_radius_doi_percentage,doi,energy_budget,epsilon,scenario,algorithm,area_length,area_width,energy_cons_fly,energy_cons_hover,wireless_technology,iterations,tours_number_avg,tours_number_std,uncovered_sensors_avg,uncovered_sensors_std,lost_data_avg,lost_data_std,seed,running_time_avg,running_time_std"
                << endl;
        file << par.num_sensors << ","
             << par.num_depots << ","
             << par.sensor_radius << ","
             << par.sensor_radius_doi_percentage << ","
             << par.doi << ","
             << par.energy_budget / 1.e+6 << "," // in MJ
             << par.epsilon << ","
             << par.scenario << ","
             << par.algorithm << ","
             << par.area_length << ","
             << par.area_width << ","
             << par.energy_cons_fly << ","
             << par.energy_cons_hover << ","
             << par.wireless_technology << ","
             << par.iterations << ","
             << tours_number_avg << ","
             << tours_number_std << ","
             << uncovered_sensors_avg << ","
             << uncovered_sensors_std << ","
             << lost_data_avg << ","
             << lost_data_std << ","
             << par.seed << ","
             << running_time_avg << ","
             << running_time_std << endl;
        file.close();
        cout << "Output saved to: " << filename << endl;
    } else {
        cerr << "Error: Unable to open file " << filename << endl;
    }
}