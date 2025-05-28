#include <matplot/matplot.h>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>



using namespace std;
using namespace matplot;
namespace py = pybind11;


vector<double> operator*(const vector<double>& vec, double scalar) {
    vector<double> result(vec.size());
    transform(vec.begin(), vec.end(), result.begin(), [scalar](double v) { return v * scalar; });
    return result;
}

vector<double> operator*(double scalar, const vector<double>& vec) {
    return vec * scalar;
}

vector<double> operator+(const vector<double>& a, const vector<double>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Wektory musza miec te sama liczbe probek!");
    }

    vector<double> result(a.size());
    transform(a.begin(), a.end(), b.begin(), result.begin(), [](double x, double y) { return x + y; });
    return result;
}

void visualize_signal(const vector<double>& sygnal, double probkowanie, double czas_od, double czas_do) {
    int liczba_probek = (czas_do - czas_od) * probkowanie;
    vector<double> czas(liczba_probek);

    // generowanie wektora czasu
    for (size_t i = 0; i < liczba_probek; ++i) {
        czas[i] = czas_od + i / probkowanie;
    }

    plot(czas, sygnal);
    xlabel("Czas (s)");
    ylabel("Wartosci Y");
    grid(true);
    show();
}

void visualize_signal_dft(const vector<complex<double>>& sygnal, double probkowanie) {
    int n = sygnal.size();
    vector<double> freq(n / 2);
    vector<double> amplitudes(n / 2);

    // Generowanie osi cz�stotliwo�ci
    for (int i = 0; i < n / 2; ++i) {
        freq[i] = i * probkowanie / n;
        amplitudes[i] = abs(sygnal[i]);
    }

    plot(freq, amplitudes);
    xlabel("Cz�stotliwo�� (Hz)");
    ylabel("Amplituda");
    title("Widmo amplitudowe (DFT)");
    grid(true);
    show();
}

vector<complex<double>> dft(const vector<double>& sygnal) {
    int n = sygnal.size();
    vector <complex<double>> wynik(n);
    for (int i = 0; i < n; ++i) {
        complex<double> suma(0, 0);
        for (int j = 0; j < n; ++j) {
            double angle = -2 * 3.14 * i * j / n;
            suma += sygnal[j] * complex<double>(cos(angle), sin(angle));
        }
        wynik[i] = suma;
    }
    return wynik;
}

vector<double> idft(const vector<complex<double>>& sygnal) {
    int n = sygnal.size();
    vector<double> wynik(n);
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0, 0);
        for (int j = 0; j < n; ++j) {
            double angle = 2 * 3.14 * j * i / n;
            sum += sygnal[j] * complex<double>(cos(angle), sin(angle));
        }
        wynik[i] = real(sum) / n;
    }
    return wynik;
}

vector<double> filtracja_1D(const vector<double>& sygnal, int k) {
    int n = sygnal.size();
    vector<double> wynik(n);
    int half = k / 2;

    for (int i = 0; i < sygnal.size(); ++i) {
        double suma = 0;
        int c = 0;
        for (int j = -half; j <= half; ++j) {
            int idx = i + j;
            if (idx >= 0 && idx < sygnal.size()) {
                suma += sygnal[idx];
                ++c;
            }
        }
        wynik[i] = suma / c;
    }
    
    return wynik;
}


py::array_t<double> filtracja_2D(const py::array_t<double>& input, int window) {
    auto buf = input.unchecked<2>();
    int rows = buf.shape(0);
    int cols = buf.shape(1);

    vector<vector<double>> temp(rows, vector<double>(cols));

    // Filtracja wierszy
    for (int r = 0; r < rows; ++r) {
        vector<double> wiersz(cols);
        for (int c = 0; c < cols; ++c)
            wiersz[c] = buf(r, c);

        auto przefiltrowany_wiersz = filtracja_1D(wiersz, window);
        for (int c = 0; c < cols; ++c)
            temp[r][c] = przefiltrowany_wiersz[c];
    }

    vector<vector<double>> wynik(rows, vector<double>(cols));

    // Filtracja kolumn
    for (int c = 0; c < cols; ++c) {
        vector<double> kolumna(rows);
        for (int r = 0; r < rows; ++r)
            kolumna[r] = temp[r][c];

        auto przefiltrowana_kolumna = filtracja_1D(kolumna, window);
        for (int r = 0; r < rows; ++r)
            wynik[r][c] = przefiltrowana_kolumna[r];
    }

    // Konwersja z vector na array
    return py::cast(wynik);
}


vector<double> generate_sin(double czestotliwosc, double probkowanie, double czas_od, double czas_do) {
    int liczba_probek = (czas_do - czas_od) * probkowanie;
    vector<double> sygnal(liczba_probek);
    for (int i = 0; i < liczba_probek; ++i) {
        double t = czas_od + i / probkowanie;
        sygnal[i] = sin(2 * 3.14 * czestotliwosc * t);
    }
    return sygnal;
}

vector<double> generate_cos(double czestotliwosc, double probkowanie, double czas_od, double czas_do) {
    int liczba_probek = (czas_do - czas_od) * probkowanie;
    vector<double> sygnal(liczba_probek);
    for (int i = 0; i < liczba_probek; ++i) {
        double t = czas_od + i / probkowanie;
        sygnal[i] = cos(2 * 3.14 * czestotliwosc * t);
    }
    return sygnal;
}

vector<double> generate_square(double czestotliwosc, double probkowanie, double czas_od, double czas_do, double amplituda, double wypelnienie) {
    int liczba_probek = (czas_do - czas_od) * probkowanie;
    vector<double> sygnal(liczba_probek);
    for (int i = 0; i < liczba_probek; ++i) {
        double t = czas_od + i / probkowanie;
        double okres = 1.0 / czestotliwosc;
        double czas_w_cyklu = fmod(t, okres);

        if (czas_w_cyklu < wypelnienie * okres) {
            sygnal[i] = amplituda;
        }
        else {
            sygnal[i] = 0;
        }
    }
    return sygnal;
}

vector<double> generate_sawtoothe(double czestotliwosc, double probkowanie, double czas_od, double czas_do, double amplituda, double offset) {
    int liczba_probek = (czas_do - czas_od) * probkowanie;
    vector<double> sygnal(liczba_probek);
    double okres = 1.0 / czestotliwosc;

    for (size_t i = 0; i < liczba_probek; ++i) {
        double t = czas_od + i / probkowanie;
        double czas_w_cyklu = fmod(t, okres);
        sygnal[i] = amplituda * (2.0 * (czas_w_cyklu / okres) - 1.0) + offset;
    }

    return sygnal;
}

vector<double> progowanie(const vector<double>& sygnal, double prog) {
    int n = sygnal.size();
    vector<double> wynik(n);
    for (int i = 0; i < n; ++i) {
        wynik[i] = ( sygnal[i] >= prog ? 1 : 0 );
    }
    return wynik;
}
PYBIND11_MODULE(signal_processing, m) {
    m.doc() = "Moduł do przetwarzania sygnałów (1D i 2D) oraz DFT";

    m.def("generate_sin", &generate_sin, "Generuje sygnał sinusoidalny");
    m.def("generate_cos", &generate_cos, "Generuje sygnał cosinusoidalny");
    m.def("generate_square", &generate_square, "Generuje sygnał prostokątny", 
          py::arg("czestotliwosc"), py::arg("probkowanie"), py::arg("czas_od"), py::arg("czas_do"), py::arg("amplituda"), py::arg("wypelnienie"));
    m.def("generate_sawtoothe", &generate_sawtoothe, "Generuje sygnał piłokształtny",
          py::arg("czestotliwosc"), py::arg("probkowanie"), py::arg("czas_od"), py::arg("czas_do"), py::arg("amplituda"), py::arg("offset"));
    m.def("progowanie", &progowanie, "Prógowanie sygnału");

    m.def("visualize_signal", &visualize_signal, "Wizualizuje sygnał czasowy");
    m.def("visualize_signal_dft", &visualize_signal_dft, "Wizualizuje widmo amplitudowe (DFT)");

    m.def("dft", &dft, "Oblicza dyskretną transformatę Fouriera (DFT)");
    m.def("idft", &idft, "Oblicza odwrotną DFT (IDFT)");

    m.def("filtracja_1D", &filtracja_1D, "Filtracja 1D średnią kroczącą");
    m.def("filtracja_2D", &filtracja_2D, "Filtracja 2D średnią kroczącą");
}