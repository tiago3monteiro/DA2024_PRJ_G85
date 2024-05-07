//
// Created by tiagomonteiro on 5/7/24.
//

#ifndef DA2024_PRJ2_G85__APPLICATION_H
#define DA2024_PRJ2_G85__APPLICATION_H


#include "Graph.h"

class Application {

public:
    Application(int i);
    void tspBacktracking();

private:
    Graph graph;
    std::vector<std::vector<float>> distanceMatrix;
    void tsp(int currPos, int n, int count, float cost, float &ans);
};


#endif //DA2024_PRJ2_G85__APPLICATION_H
