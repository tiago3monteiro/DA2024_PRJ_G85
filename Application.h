#ifndef DA2024_PRJ2_G85__APPLICATION_H
#define DA2024_PRJ2_G85__APPLICATION_H


#include <set>
#include "Graph.h"

class Application {

public:
    Application(int i);
    /**
     * @brief Determines an optimal solution to the Travelling Salesman Problem
     *
     * This function employs a brute-force approach by testing all possible hamiltonian cycles of the graph
     * and displaying the minimum cost found. It keeps track of the current minimum cost and updates it
     * as it finds lower costs. After testing a particular path, it backtracks and starts testing another path.
     *
     * **Complexity**: O(V!), where V is the number of vertices in the graph
     *
      */
    void tspBacktracking();
    /**
     * @brief Polynomial approximation algorithm that relies on the triangular inequality principle
     *
     * This function uses a polynomial-time triangular approximation to solve the Travelling Salesman Problem.
     * It constructs a Minimum Spanning Tree (MST) using Prim's algorithm, performs a preorder traversal of the
     * MST to generate a TSP tour, and calculates the total cost. If the graph lacks geographical coordinates,
     * the function is not applicable.
     *
     * **Complexity**: O(V + E), where V is the number of vertices and E is the number of edges in the graph
     *
      */
    void tspTriangular();
    /**
     * @brief Implements the nearest neighbor heuristic to solve the Travelling Salesman Problem
     *
     * Starting from the first vertex, it iteratively selects the nearest unvisited vertex until all vertices are
     * visited, then returns to the starting vertex to complete the tour. If the graph lacks geographical
     * coordinates, the function is not applicable.

     *
     * **Complexity**: O(V^2), where V is the number of vertices in the graph
     *
      */
    void tspNearestNeighbor();
    /**
     * @brief Uses Christofides' algorithm to solve the Travelling Salesman Problem.
     *
     * It constructs a Minimum Spanning Tree (MST), finds a minimum weight perfect matching for odd-degree vertices,
     * and combines these to form an Eulerian circuit. The Eulerian circuit is then converted into a Hamiltonian circuit
     * (TSP tour), and the total cost is calculated. If the graph lacks geographical coordinates, the function is
     * not applicable.

     *
     * **Complexity**: O(E * V^2), where V is the number of vertices and E is the number of edges in the graph
     *
      */
    void tspChristofides();
     /**
     * @brief Uses a nearest neighbor approach to solve the Travelling Salesman Problem in real-world scenarios.
     *
     * Starting from a given source vertex, it iteratively selects the nearest unvisited vertex
     * until all vertices are visited, then returns to the source to complete the tour. The function calculates
     * and displays the total distance of the tour.

     *
     * **Complexity**: O(V^2), where V is the number of vertices in the graph
     *
      */
    void tspRealWorld(int source);

private:
    Graph graph;
    std::vector<std::vector<float>> distanceMatrix;
    std::vector<std::vector<bool>> visited; // Matrix to track visited edges
    std::vector<std::vector<int>> mst; // To store the MST as adjacency list
    std::vector<int> tspTour; // To store the TSP tour sequence
    int graphType;

    /**
     * @brief Auxiliary function to backtracking algorithm
     *
     * It recursively explores all possible Hamiltonian cycles, updating the minimum cost and path whenever a
     * lower-cost tour is found. The function backtracks by unmarking vertices and removing them from the
     * path to explore alternative routes.

     *
     * **Complexity**: O(V!), where V is the number of vertices in the graph
     *
      */
    void tspBacktrackingAux(int currPos, int n, int count, float cost, float &ans, std::vector<int>& path); //auxiliar function for backtracking
    /**
     * @brief Prim's algorithm to find Minimum Spanning Tree (MST) of graph
     *
     * Uses a priority queue to select the minimum weight edge at each step, ensuring each vertex is included in the MST.
     * The function calculates the total weight of the MST and updates the MST structure.

     *
     * **Complexity**: O(E * log(V)), where V is the number of vertices and E the number of edges in the graph
     *
      */
    void primMST();
     /**
     * @brief Performs a preorder traversal starting from a given root node
     *
     * It recursively visits each node and adds it to the TSP tour if it hasn't been visited before.
     * The traversal ensures that every vertex in the MST is included in the tour.

     *
     * **Complexity**: O(V), where V is the number of vertices in the graph
     *
      */
    void preorderTraversal(int root , std::vector<bool> &visited);
    /**
     * @brief Determines Haversine distance between two points
     *
     * It makes use of the geographical coordinates of the two vertices and applies the formula to calculate
     * the distance between the two.

     *
     * **Complexity**: O(log(n))
     *
      */
    double haversineDistance(Vertex *v1, Vertex *v2);
    /**
     * @brief Implements the Blossom Algorithm to find perfect matching
     *
     * It iteratively constructs augmenting paths by using a priority queue to find the shortest path from an
     * unmatched vertex to any other vertex. The algorithm continues until all odd vertices are matched.
     * The matching edges are then extracted and returned.

     *
     * **Complexity**: O(E * V^2), where V is the number of vertices and E the number of edges in the graph
     *
      */
    std::vector<std::pair<int, int>> blossomAlgorithm(const std::vector<int>& oddVertices);
    /**
     * @brief Finds an Eulerian circuit in a graph
     *
     * It traverses each edge of the MST, marking it as visited to ensure that no edge is traversed twice.
     * The function then adds the vertices to the circuit vector in a depth-first manner.

     *
     * **Complexity**: O(E), where E is the number of edges in the graph
     *
      */
    void findEulerianCircuit(int u, std::vector<int> &circuit);
    void resetGraph();

};


#endif //DA2024_PRJ2_G85__APPLICATION_H
