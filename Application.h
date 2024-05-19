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
    int graphType;

    void tspBacktrackingAux(int currPos, int n, int count, float cost, float &ans, std::vector<int>& path); //auxiliar function for backtracking
    void primMST();
    void preorderTraversal(int root , std::vector<bool> &visited);
    double haversineDistance(Vertex *v1, Vertex *v2);
    std::vector<std::pair<int, int>> blossomAlgorithm(const std::vector<int>& oddVertices);
    void findEulerianCircuit(int u, std::vector<int> &circuit);
    void resetGraph();

};


#endif //DA2024_PRJ2_G85__APPLICATION_H
