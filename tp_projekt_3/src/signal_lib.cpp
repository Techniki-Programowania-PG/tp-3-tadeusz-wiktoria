#include "signal_lib.h"
#include <cmath>
#include "/matplotplusplus/source/matplot/matplot.h"
#include <iostream>
#include <filesystem>
#include "signal_lib.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <complex>

using namespace std;
using namespace matplot;

std::vector<double> add_noise(const std::vector<double>& signal, double amp) {
    std::vector<double> noisy(signal);
    std::default_random_engine gen;
    std::normal_distribution<double> dist(0.0, amp);
    for (auto& val : noisy) {
        val += dist(gen);
    }
    return noisy;
}

std::vector<int> thresholding(const std::vector<double>& signal, double threshold) {
    std::vector<int> result(signal.size());
    for (size_t i = 0; i < signal.size(); ++i)
        result[i] = signal[i] > threshold ? 1 : 0;
    return result;
}

std::vector<double> low_pass_filter(const std::vector<double>& signal, double cutoff) {
    int kernel_size = static_cast<int>(signal.size() * cutoff);
    kernel_size = std::max(1, kernel_size);
    std::vector<double> result(signal.size(), 0.0);
    for (size_t i = 0; i < signal.size(); ++i) {
        for (int j = 0; j < kernel_size; ++j) {
            int idx = i + j - kernel_size / 2;
            if (idx >= 0 && idx < signal.size()) {
                result[i] += signal[idx];
            }
        }
        result[i] /= kernel_size;
    }
    return result;
}

std::vector<double> derivative(const std::vector<double>& signal) {
    std::vector<double> result(signal.size());
    result[0] = 0;
    for (size_t i = 1; i < signal.size(); ++i)
        result[i] = signal[i] - signal[i - 1];
    return result;
}

std::vector<double> edge_detect(const std::vector<double>& signal) {
    std::vector<double> result(signal.size(), 0.0);
    for (size_t i = 1; i + 1 < signal.size(); ++i)
        result[i] = signal[i + 1] - signal[i - 1];
    return result;
}

std::vector<double> remove_low_freq(const std::vector<double>& signal, double ratio) {
    size_t N = signal.size();
    std::vector<std::complex<double>> spectrum(N);
    for (size_t k = 0; k < N; ++k) {
        for (size_t n = 0; n < N; ++n) {
            double angle = -2 * M_PI * k * n / N;
            spectrum[k] += signal[n] * std::complex<double>(cos(angle), sin(angle));
        }
    }
    int k = static_cast<int>(ratio * N);
    for (int i = 0; i < k; ++i) spectrum[i] = 0;
    for (int i = N - k; i < N; ++i) spectrum[i] = 0;

    std::vector<double> result(N);
    for (size_t n = 0; n < N; ++n) {
        std::complex<double> sum = 0;
        for (size_t k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            sum += spectrum[k] * std::complex<double>(cos(angle), sin(angle));
        }
        result[n] = sum.real() / N;
    }
    return result;
}



std::string plot_signal(const std::vector<double>& signal, const std::string& filename = "my_plot.png", bool show_plot = true) {
    auto fig = figure();
    plot(signal);

    xlabel("Samples");
    ylabel("Amplitude");
    title("Signal Plot");

    fig->save(filename);  // Сохраняем файл
    std::string full_path = std::filesystem::absolute(filename).string();

    if (show_plot) {
        show();
    }

    return full_path;
}

std::vector<double> generate_sin(double frequency, double sample_rate, int samples) {
    std::vector<double> signal(samples);
    for (int i = 0; i < samples; ++i)
        signal[i] = std::sin(2 * 3.14 * frequency * i / sample_rate);
    return signal;
}

std::vector<double> generate_cos(double frequency, double sample_rate, int samples) {
    std::vector<double> signal(samples);
    for (int i = 0; i < samples; ++i)
        signal[i] = std::cos(2 * 3.14 * frequency * i / sample_rate);
    return signal;
}

std::vector<double> generate_square(double frequency, double sample_rate, int samples) {
    std::vector<double> signal(samples);
    for (int i = 0; i < samples; ++i) {
        double t = i / sample_rate;
        signal[i] = std::sin(2 * 3.14 * frequency * t) >= 0 ? 1.0 : -1.0;
    }
    return signal;
}

std::vector<double> generate_sawtooth(double frequency, double sample_rate, int samples) {
    std::vector<double> signal(samples);
    for (int i = 0; i < samples; ++i) {
        double t = i / sample_rate;
        signal[i] = 2.0 * (t * frequency - floor(t * frequency + 0.5));
    }
    return signal;
}
