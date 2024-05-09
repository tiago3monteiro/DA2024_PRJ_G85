#include <iostream>
#include "Application.h"

int main() {
    Application app(8);
    std::cout << "parse done";
    std::cout.flush();
    //app.tspBacktracking();
    app.tspTriangular();
    return 0;
}
