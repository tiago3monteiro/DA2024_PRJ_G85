//
// Created by tiagomonteiro on 5/7/24.
//

#ifndef DA2024_PRJ2_G85__APPLICATION_H
#define DA2024_PRJ2_G85__APPLICATION_H


#include <set>
#include "Graph.h"

class Application {

public:
    Application(int i);
    void tspBacktracking();
    void tspTriangular();

private:
    Graph graph;
    std::vector<std::vector<float>> distanceMatrix;
    std::vector<std::vector<int>> mst; // To store the MST as adjacency list
    std::vector<int> tspTour; // To store the TSP tour sequence

    void tspBacktrackingAux(int currPos, int n, int count, float cost, float &ans); //auxiliar function for backtracking
    void primMST();
    void preorderTraversal(int root , std::vector<bool> &visited);
    int minKey(const std::vector<float> &key, const std::vector<bool> &mstSet);
    double haversineDistance(Vertex *v1, Vertex *v2);
};


#endif //DA2024_PRJ2_G85__APPLICATION_H
