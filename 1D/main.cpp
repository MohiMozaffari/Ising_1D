#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;


/**
 * @brief Function to write vectors to a CSV file
 * 
 * @param temp temprature vector
 * @param Mag magnetization vector
 * @param Energy energy vector
 * @param filename file name
 */
void writeVectorsToCSV(const vector<double>& temp, const vector<double>& Mag, const vector<double>& Energy, const string& filename) {
    ofstream file(filename);

    if (file.is_open()) {
        // Write column headers
        file << "Temperature,Magnetization,Energy\n";

        // Write data to the file
        size_t size = temp.size();
        for (size_t i = 0; i < size; ++i) {
            file << temp[i] << "," << Mag[i] << "," << Energy[i] << "\n";
        }

        // Close the file
        file.close();
        cout << "File saved successfully!" << endl;
    }
    else {
        cout << "Error opening the file!" << endl;
    }
}


int main()
{
    srand(time(NULL));
    const int N = 64;
    const int nens = 50;
    const double minTemp = 0.0001;
    const double maxTemp = 4.0;
    const int numTemps = 40;

    vector<double> temp(numTemps);
    double step = (maxTemp - minTemp) / (numTemps - 1);
    for (int i = 0; i < numTemps; i++) {
        temp[i] = minTemp + step * i;
    }

    vector<double> Energy(numTemps);
    vector<double> Mag(numTemps);

    auto start = high_resolution_clock::now();
    vector<int> lattice (N, 0);
    for (int i = 0; i < N; i++)
    {
        lattice[i] = (rand()%2 ? 1: -1);
        //lattice[i] = 1;
    }
    for (int t = 0; t < numTemps; t++) {
        //static const vector<double> expTable = {1.0, exp(-4.0 / temp[t]), exp(-8.0 / temp[t])};
        double avgE = 0.0;
        double avgM = 0.0;
        for (int ens = 0; ens < nens; ens++) {

            for (int i = 0; i < pow(N,3); i++) {

                int n = (rand() % (N-1));
                int dE = -2*lattice[n] * (lattice[(n + 1) % N] + lattice[(n - 1 + N) % N]);
                dE -= lattice[n] * (lattice[(n + 2) % N] + lattice[(n - 2 + N) % N]);
                if (dE <= 0 || ((rand() / RAND_MAX) < exp(-dE / temp[t]))) {
                    lattice[n] *= -1;
                }

            }

            double mag = 0.0;
            for (const auto spin : lattice) {
                mag += spin;
            }

            double energy = 0.0;
            for (size_t n = 0; n < lattice.size(); n++) {
                    energy -= lattice[n]*(lattice[(n+1) % lattice.size()] + lattice[(n+2) % lattice.size()]);
            }

            avgE += energy/N;
            avgM += abs(mag)/N;
        }

        avgE /= nens;
        avgM /= nens;

        Energy[t] = avgE;
        Mag[t] = avgM;
    }
    
    string filePrefix = "C:\\Users\\Asus\\Documents\\Programming\\Ising_1D\\first_neighbor\\Ising_1DII_";
    string fileExtension = ".csv";
    string fileName = filePrefix+ to_string(N) + fileExtension;

    writeVectorsToCSV(temp, Mag, Energy, fileName);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Simulation completed in " << duration.count() << " seconds." << endl;

    return 0;
}