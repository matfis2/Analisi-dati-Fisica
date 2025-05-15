#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

int main() {
    string inputFileName, outputFileName;
    double min,max;
    double pr,vo,te;
    vector<double> p,v,t;

    // Richiesta all'utente dei nomi dei file e della soglia
    cout << "Nome del file di input: ";
    getline(cin, inputFileName);
    cout << "Nome del file di output: ";
    getline(cin, outputFileName);
    cout << "Inserisci la soglia minima: ";
    cin >> min;
    cout<<"inserisci la soglia massima di volume:";
    cin>>max;

    ifstream inputFile(inputFileName);
    if (!inputFile) {
        cerr << "Errore nell'apertura del file di input.\n";
        return 1;
    }

    while(inputFile>>pr>>vo>>te){
        v.push_back(vo);
        p.push_back(pr);
        t.push_back(te);
    }
    ofstream outputFile(outputFileName);
    if (!outputFile) {
        cerr << "Errore nella creazione del file di output.\n";
        return 1;
    }
    for(int i=0;i<v.size();i++){
        if(v.at(i)>min && v.at(i)<max){
            outputFile << p.at(i) << ' ' << v.at(i) << ' ' << t.at(i) << '\n';
        }
    }
    cout << "Filtraggio completato. Dati salvati in '" << outputFileName << "'.\n";
    inputFile.close();
    outputFile.close();
    return 0;
}