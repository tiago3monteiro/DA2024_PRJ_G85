#include <iostream>
#include "Application.h"

int main() {
    Application app(8);
    std::cout << "parse done"<<std::endl;
    std::cout.flush();
    //app.tspBacktracking();
    app.tspTriangular();
    app.tspNearestNeighbor();
    app.tspChirstofides();
    return 0;
}
