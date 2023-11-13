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
 * @brief Function to initialize the random 1D lattice btw -1 and 1
 * 
 * @param N number of lattice
 * @return vector<int> lattice 1d 
 */
vector<int> initializer(int N)
{
    vector<int> lattice (N, 0);
    for (int i = 0; i < N; i++)
    {
        lattice[i] = (rand()%2 ? 1: -1);
    }
    
    return lattice;
}

/**
 * @brief Function to calculate the Hamiltonian at a given position
 * 
 * @param n position of spin
 * @param lattice list of spin
 * @return double energy
 */
double hamiltonian(int n, const vector<int>& lattice) {
    int N = lattice.size();
    double E = lattice[n] * (lattice[(n + 1) % N] + lattice[(n - 1 + N) % N]);
    E += lattice[n] * (lattice[(n + 2) % N] + lattice[(n - 2 + N) % N]);
    return E;
}

/**
 * @brief Function to flip the spin at a given position if the energy change is favorable
 * 
 * @param dE energy diff
 * @param n position of spin
 * @param lattice list of spin
 * @param ex boltzmann probability
 */
void check_flip(int n, double dE, vector<int>& lattice, double ex) {

    if (dE <= 0 || ((rand() / RAND_MAX) < ex)) {
        lattice[n] *= -1;
    }
}


/**
 * @brief Function to calculate the magnetization of the lattice
 * 
 * @param lattice list of spin
 * @return double magnetization
 */
double magnetization(const vector<int> &lattice) {
    double mag = 0.0;
    for (const auto spin : lattice) {
        mag += spin;
    }
    return (abs(mag)/lattice.size());
}


/**
 * @brief Function to calculate the energy of the lattice
 * 
 * @param lattice list of spin
 * @return double energy
 */
double energy(const vector<int>& lattice) {
    double energy = 0.0;
    for (size_t n = 0; n < lattice.size(); n++) {
        energy += hamiltonian(n, lattice);
    }
    return (energy / lattice.size()) ;
}


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
    const double minTemp = 0.01;
    const double maxTemp = 20.0;
    const int numTemps = 40;

    vector<double> temp(numTemps);
    double step = (maxTemp - minTemp) / (numTemps - 1);
    for (int i = 0; i < numTemps; i++) {
        temp[i] = minTemp + step * i;
    }

    vector<double> Energy(numTemps);
    vector<double> Mag(numTemps);

    auto start = high_resolution_clock::now();

 for (int t = 0; t < numTemps; t++) {
        static const vector<double> expTable = {1.0, exp(-4.0 / temp[t]), exp(-8.0 / temp[t]), exp(-12.0 / temp[t])};
        double avgE = 0.0;
        double avgM = 0.0;

        for (int ens = 0; ens < nens; ens++) {
            vector<int> lattice = initializer(N);

            for (int i = 0; i < pow(N,3); i++) {

                int pos = (rand() % N);
                double dE = -2 * hamiltonian(pos, lattice);
                check_flip(dE, pos, lattice, expTable[static_cast<int>(-dE) / 4]);
            }

            double E = energy(lattice);
            double M = magnetization(lattice);

            avgE += E;
            avgM += M;
        }

        avgE /= nens;
        avgM /= nens;

        Energy[t] = avgE;
        Mag[t] = avgM;
    }
    
    string filePrefix = "C:\\Users\\Asus\\Documents\\Programming\\Ising_1D\\second_neighbor\\Ising_1D_";
    string fileExtension = ".csv";
    string fileName = filePrefix+ to_string(N) + fileExtension;

    writeVectorsToCSV(temp, Mag, Energy, fileName);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Simulation completed in " << duration.count() << " seconds." << endl;

    return 0;
}




