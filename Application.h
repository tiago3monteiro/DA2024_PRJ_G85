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
    void tspNearestNeighbor();
    void tspChristofides();
    void tspRealWorld(int source);

private:
    Graph graph;
    std::vector<std::vector<float>> distanceMatrix;
    std::vector<std::vector<bool>> visited; // Matrix to track visited edges
    std::vector<std::vector<int>> mst; // To store the MST as adjacency list
    std::vector<int> tspTour; // To store the TSP tour sequence

    void tspBacktrackingAux(int currPos, int n, int count, float cost, float &ans); //auxiliar function for backtracking
    void primMST();
    void preorderTraversal(int root , std::vector<bool> &visited);
    int minKey(const std::vector<float> &key, const std::vector<bool> &mstSet);
    double haversineDistance(Vertex *v1, Vertex *v2);
    void findEulerianCircuit(int u, std::vector<int> &circuit);

    bool hamiltonianUtil(int v, std::vector<int> &path, std::vector<bool> &visited, int &count);

    bool isHamiltonian(const std::vector<int> &path);
};


#endif //DA2024_PRJ2_G85__APPLICATION_H
