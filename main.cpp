#include <iostream>
#include "Menu.h"

int main() {
    Menu menu;
    while (menu.run());
    std::cout << "PROGRAM ENDED";
    return 0;
}
