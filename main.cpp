#include <iostream>
#include "Application.h"

int main() {
    Application app(4);
    std::cout << "parse done";
    std::cout.flush();
   // app.tspBacktracking();
    app.tspTriangular();
    return 0;
}
