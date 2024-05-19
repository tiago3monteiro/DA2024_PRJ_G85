#ifndef DA2024_PRJ2_G85__MENU_H
#define DA2024_PRJ2_G85__MENU_H

#include "Application.h"

class Menu {
public:
    void clearTerminal();
    bool run();
    int processInput(const std::string input);
private:
    int key = -1;
    int mainMenu();


    bool tspMenu(Application app);
};


#endif //DA2024_PRJ2_G85__MENU_H
