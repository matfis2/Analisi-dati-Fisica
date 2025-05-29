#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::string nomeFile;
    std::cout << "Inserisci il nome del file: ";
    std::getline(std::cin, nomeFile);  // Per gestire anche nomi con spazi

    // Legge tutto il contenuto del file in una stringa
    std::ifstream fileInput(nomeFile);
    if (!fileInput) {
        std::cerr << "Errore nell'apertura del file in lettura.\n";
        return 1;
    }

    std::string contenuto((std::istreambuf_iterator<char>(fileInput)),
                          std::istreambuf_iterator<char>());
    fileInput.close();

    // Sostituisce le virgole con i punti
    for (char& c : contenuto) {
        if (c == ',') {
            c = '.';
        }
    }

    // Sovrascrive il file con il nuovo contenuto
    std::ofstream fileOutput(nomeFile);
    if (!fileOutput) {
        std::cerr << "Errore nell'apertura del file in scrittura.\n";
        return 1;
    }

    fileOutput << contenuto;
    fileOutput.close();

    std::cout << "Sostituzione completata con successo.\n";
    return 0;
}
