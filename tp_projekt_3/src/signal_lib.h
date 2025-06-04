#pragma once
#include <vector>
#include <string>

std::vector<double> generate_sin(double frequency, double sample_rate, int samples);
std::vector<double> generate_cos(double frequency, double sample_rate, int samples);
std::vector<double> generate_square(double frequency, double sample_rate, int samples);
std::vector<double> generate_sawtooth(double frequency, double sample_rate, int samples);
std::vector<double> add_noise_to_signal(const std::vector<double>& signal, double amp = 0.5);
std::vector<int> thresholding_to_signal(const std::vector<double>& signal, double threshold);
std::vector<double> low_pass_filter_to_signal(const std::vector<double>& signal, double cutoff = 0.05);
std::vector<double> derivative_to_signal(const std::vector<double>& signal);
std::vector<double> edge_detect_to_signal(const std::vector<double>& signal);
std::vector<double> remove_low_freq_to_signal(const std::vector<double>& signal, double ratio = 0.1);

std::string plot_signal(const std::vector<double>& signal);
