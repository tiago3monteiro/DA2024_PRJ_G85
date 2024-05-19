#include <limits>
#include "Menu.h"

bool Menu::run() {
    int key = mainMenu();
    if(!key) return false;
    if(key>8) return false;
    Application app(key);
    while(tspMenu(app));
    return true;
}

int Menu::processInput(const std::string input) {
    try {
        int option = std::stoi(input);
        if (option < 0 || input.size() > 2) {
            throw std::invalid_argument("");
        }
        return option;
    }
    catch (const std::invalid_argument& ia) {
        key = -1;
        std::cout << "\nError, please input a valid option." << std::endl;
        return -1;
    }
}

int Menu::mainMenu() {
    std::cout << "================TSP Choose Graph to be processed================" << std::endl;
    std::cout << "1. Small Graph - tourism" << std::endl;
    std::cout << "2. Small Graph - stadiums" << std::endl;
    std::cout << "3. Small Graph - shipping" << std::endl;
    std::cout << "4. Large Graph 1" << std::endl;
    std::cout << "5. Large Graph 2" << std::endl;
    std::cout << "6. Large Graph 3" << std::endl;
    std::cout << "7. Medium Graph 25" << std::endl;
    std::cout << "8. Medium Graph  500" << std::endl;
    std::cout << "0. Exit" << std::endl;

    std::string in;
    std::cout << "> ";
    std::cin >> in;
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    key = processInput(in);
    return key;
}

bool Menu::tspMenu(Application app) {
    std::cout << "================Choose Heuristic or Algorithm================" << std::endl;
    std::cout << "1. Backtracking" << std::endl;
    std::cout << "2. 2-approximation" << std::endl;
    std::cout << "3. Nearest Neighbor" << std::endl;
    std::cout << "4. Christofides-Serdyukov" << std::endl;
    std::cout << "5. Real World" << std::endl;
    std::cout << "6. Exit" << std::endl;
    std::string in;
    std::cout << "> ";
    std::cin >> in;
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    key = processInput(in);
    switch (key) {
        case 1:
            app.tspBacktracking();
            break;
        case 2:
            app.tspTriangular();
            break;
        case 3:
            app.tspNearestNeighbor();
            break;
        case 4:
            app.tspChristofides();
            break;
        case 5: {
            std::string source;
            std::cout << "Please enter the source for the real world TSP: ";
            std::cin >> source;
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            app.tspRealWorld(std::stoi(source));
            break;
        }
        case 6:
            return false;
        default:
            std::cout << "Not valid!" <<std::endl;
            break;
    }
    return true;
}
