//
// Created by tiagomonteiro on 5/7/24.
//

#include "Application.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <array>
#include <climits>

#include <cmath> // For math functions
Application::Application(int i) {
    std::array<std::string, 8> filenames = {
            "../graphs/tourism.csv",
            "../graphs/stadiums.csv",
            "../graphs/shipping.csv",
            "../graphs/graph1/",
            "../graphs/graph2/",
            "../graphs/graph3/",
            "../graphs/graph4-25/",
            "../graphs/graph5-500/"

    };
    //PARSE TOY GRAPHS:
    if(i > 0 && i<4)
    {
        std::string filepath =  filenames[i-1];
        std::ifstream in(filepath);
        std::string line;
        std::getline(in,line,'\n'); //remove first line

        while(std::getline(in,line,'\n')) {
            //line = line.substr(0, line.length() - 1);
            std::istringstream iss(line);
            std::string word;
            std::vector<std::string> saved;
            while (std::getline(iss, word, ',')) saved.push_back(word);
            Vertex *sourceVertex = this->graph.findVertex(stoi(saved[0])); //source
            Vertex *destVertex = this->graph.findVertex(stoi(saved[1])); //destiny

            if (sourceVertex == nullptr) { //vertex does not exist yet
                sourceVertex = new Vertex(stoi(saved[0]));
                graph.addVertex(sourceVertex);
            }
            if (destVertex == nullptr) { //vertex does not exist yet
                destVertex = new Vertex(stoi(saved[1]));
                graph.addVertex(destVertex);
            }
            graph.addEdge(stoi(saved[0]), stoi(saved[1]), stof(saved[2]));
            graph.addEdge(stoi(saved[1]), stoi(saved[0]), stof(saved[2]));//add reverse

        }
    }

    else {
        std::string edges =  filenames[i-1] + "/edges.csv";
        std::string nodes =  filenames[i-1] + "/nodes.csv";
        std::ifstream in(nodes);
        std::string line;
        std::getline(in,line,'\n'); //remove first line
        //parse nodes:

        while(std::getline(in,line,'\n')) {

            std::istringstream iss(line);
            std::string word;
            std::vector<std::string> saved;
            while (std::getline(iss, word, ',')) saved.push_back(word);
            auto vertex = new Vertex(stoi(saved[0]));
            vertex->setLon(stod(saved[1]));
            vertex->setLat(stod(saved[2]));
            graph.addVertex(vertex);
        }
        int k = 0;
        //parse edges:
        std::ifstream in1(edges);
        std::getline(in1,line,'\n'); //remove first line
        while(std::getline(in1,line,'\n')) {
            std::cout << "line " << k << std::endl;
            k++;

            std::istringstream iss1(line);
            std::string word;
            std::vector<std::string> saved;
            while (std::getline(iss1, word, ',')) saved.push_back(word);
            graph.addEdge(stoi(saved[0]),stoi(saved[1]),stof(saved[2]));
            graph.addEdge(stoi(saved[1]),stoi(saved[0]),stof(saved[2]));

        }
    }

    //distance matrix:
    int n = graph.getVertexSet().size();
    visited.resize(n, std::vector<bool>(n, false));
    distanceMatrix.assign(n, std::vector<float>(n, std::numeric_limits<float>::infinity()));

    for(auto vertex: graph.getVertexSet()) {
        distanceMatrix[vertex->getCode()][vertex->getCode()] = 0;
    }
    for(auto vertex: graph.getVertexSet()) {
        for (auto edge: vertex->getAdj()) {
            distanceMatrix[edge->getOrig()->getCode()][edge->getDest()->getCode()] = edge->getCapacity();
            distanceMatrix[edge->getDest()->getCode()][edge->getOrig()->getCode()] = edge->getCapacity();//reverse

        }
    }
    graph.setDistanceMatrix(distanceMatrix);


}
//optimal solution:
//2600 -tourism
//86.7 - shipping
//341 - stadiums

//T2.1 .................................................................................................
void Application::tspBacktracking() {

    auto start = std::chrono::steady_clock::now();

    //setup:
   // n -> Number of nodes (vertices) in the graph
   // set all vertex but the source one as unvisited, ans as infinte!
    int n = distanceMatrix.size();
    for(auto vertex: graph.getVertexSet())
        vertex->setVisited(false);
    graph.findVertex(0)->setVisited(true);
    float ans = std::numeric_limits<float>::max();

    //Initiates the recursive TSP (Travelling Salesman Problem):
    // currPos -> function tsp starting from node 0,
    // n -> with the total number of nodes,
    // count -> initial count of visited nodes (1 for the starting node),
    // initial cost -> (0.0f)
    // and reference to ans.
    tspBacktrackingAux(0, n, 1, 0.0f, ans);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Lowest costing path starting at 0: " << ans << std::endl;
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;

}

void Application::tspBacktrackingAux(int currPos, int n, int count, float cost, float& ans) {

    // Base case:
    //checks if count == n, each means all vertices have been visited
    //and checks if there's a direct link between the current vertex and vertex 0 (source)
    //if it does the ans value will be updated IF ans > cost + distanceMatrix[currPos][0]
    //which means this new hamiltonian cycle has a lower cost then the one previously discovered.
    if (count == n && distanceMatrix[currPos][0] > 0) {
        ans = std::min(ans, cost + distanceMatrix[currPos][0]);
        return;
    }

    // Explore all adjacent nodes (vertices) from the current position (currPos):
    // we go through all edges of our current vertex to check its destination vertices
    //if this vertex hasn't been visited yet and exist a direct connection to the next vertex
    //we will mark it has visited and use recursion but now count -> count +1, cost -> cost + distanceMatrix[currPos][nextVertex]
   // and then backtracks by marking the destination vertex as unvisited after the recursive call.
    for (auto edge : graph.findVertex(currPos)->getAdj()) {
        int nextVertex = edge->getDest()->getCode();
        if (!edge->getDest()->isVisited() && distanceMatrix[currPos][nextVertex] > 0) {
            edge->getDest()->setVisited(true);
            tspBacktrackingAux(nextVertex, n, count + 1, cost + distanceMatrix[currPos][nextVertex], ans);
            edge->getDest()->setVisited(false); //backtrack
        }
    }
}
//grafo 1 -> 1.14207e+06
//execution time 96 milliseconds

//grafo 2 -> 529 203/528 965/2.05874e+06
//execution time -> 585803/492615/2091 milliseconds -> 8.21 minutes
//T2.2 .................................................................................................
void Application::tspTriangular() { //Complexity – O(V+E) so it is polynomial on the size of G

    auto start = std::chrono::steady_clock::now();

    primMST(); // Construct the MST using Prim's algorithm
    std::vector<bool> visited(distanceMatrix.size(), false);
    preorderTraversal(0, visited); // Perform preorder traversal of MST

    //print outcome:
    std::cout << "TSP Tour Sequence: ";
    float totalCost = 0.0f;
    for (size_t i = 0; i < tspTour.size(); i++) {
        int current = tspTour[i];
        int next = tspTour[(i + 1) % tspTour.size()]; // Wrap around to the beginning

        if (distanceMatrix[current][next]  >= INT_MAX) {
            Vertex* source = graph.findVertex(current);
            Vertex* destiny = graph.findVertex(next);
            auto distance = haversineDistance(  source, destiny);
            totalCost += distance;
        }
        else totalCost += distanceMatrix[current][next];
    }
    std::cout << "Total Cost: " << totalCost << std::endl;
    //time count:
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
}

int Application::minKey(const std::vector<float>& key, const std::vector<bool>& mstSet) {
    float min = std::numeric_limits<float>::max();
    int min_index = -1;

    for (size_t v = 0; v < key.size(); v++) {
        if (!mstSet[v] && key[v] < min) {
            min = key[v];
            min_index = v;
        }
    }

    return min_index;
}
void Application::primMST() {
    int n = distanceMatrix.size(); // Number of nodes
    std::vector<bool> mstSet(n, false); // Visited vertices
    std::vector<float> key(n, std::numeric_limits<float>::max()); // Lowest edge weight connecting vertex to the MST
    std::vector<int> parent(n, -1); // Predecessor of vertex in the MST

    key[0] = 0; // Start with vertex 0 as the root

    // Construct the MST using Prim's algorithm
    for (int count = 0; count < n - 1; count++) {
        int u = minKey(key, mstSet); // Get the closest unprocessed node
        mstSet[u] = true; // Mark vertex u as visited

        // Update key values and parent pointers for adjacent vertices
        for (int v = 0; v < n; v++) {
            if (distanceMatrix[u][v] > 0 && !mstSet[v] && distanceMatrix[u][v] < key[v]) {
                parent[v] = u; // Update parent of vertex v
                key[v] = distanceMatrix[u][v]; // Update the key value for vertex v
            }
        }
    }

    // Calculate the total weight of the MST
    float totalMSTWeight = 0.0f;
    for (int i = 1; i < n; i++) {
        totalMSTWeight += key[i]; // Sum up the key values (MST edge weights)
    }

    // Construct MST adjacency list for preorder traversal
    mst.resize(n);
    for (int v = 1; v < n; v++) {
        mst[parent[v]].push_back(v);
        mst[v].push_back(parent[v]);
    }

    // Print the total weight of the MST
    std::cout << "Total weight of MST: " << totalMSTWeight << std::endl;
}


void Application::preorderTraversal(int root, std::vector<bool>& visited) {

    if (!visited[root]) {
        tspTour.push_back(root);
        visited[root] = true;

        for (int neighbor : mst[root]) {
            if (!visited[neighbor]) {
                preorderTraversal(neighbor, visited);
            }
        }
    }
}



double Application::haversineDistance(Vertex* v1, Vertex* v2) {
    // Radius of the Earth in kilometers
    constexpr float R = 6371.0;

    // Convert latitude and longitude from degrees to radians
    double lat1 = v1->getLat() * M_PI / 180.0;
    double lon1 = v1->getLon() * M_PI / 180.0;
    double lat2 = v2->getLat() * M_PI / 180.0;
    double lon2 = v2->getLon() * M_PI / 180.0;

    // Calculate differences in latitude and longitude
    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    // Haversine formula to calculate distance
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    // Distance in kilometers
    double distance = R * c;
    return distance;
}


//5.50069e+06 units -> graph 2
void Application::tspNearestNeighbor() {

    auto start = std::chrono::steady_clock::now();
    int n = distanceMatrix.size();
    int visitedCount = 0;
    Vertex* next;
    Vertex* previous;
    std::vector<bool> visited(n, false);
    double totalDistance = 0.0;

    // Start from the first vertex
    Vertex* current = graph.findVertex(0);
    current->setVisited(true);
    visitedCount++;

    while (visitedCount < n) { //lets visit all nodes:
        bool pathExists = false;
        double minDistance = std::numeric_limits<double>::max();

        for(auto edge: current->getAdj()) { //check adjcent nodes to look for the best local choice
            if(!edge->getDest()->isVisited() && edge->getCapacity() < minDistance) {
                next = edge->getDest();
                minDistance = edge->getCapacity();
                pathExists = true; //exists at least one!
            }
        }

        if(!pathExists) {
            for(auto vertex: graph.getVertexSet()) {
                if(!vertex->isVisited()) {
                    double capacity = haversineDistance(vertex,current);
                    if(capacity < minDistance) {
                        next = vertex;
                        minDistance = capacity;
                    }
                }
            }
        }
        visitedCount++;
        previous = current;
        current = next;
        current->setVisited(true);
        totalDistance += minDistance;
    }
    // Add the distance back to the starting vertex to complete the tour
    bool path = false;
    for(auto edge: current->getAdj()) {
        if(edge->getDest()->getCode() == 0) {
            totalDistance += edge->getCapacity();
            path = true;
        }
    }
    if(!path)
    totalDistance += haversineDistance(previous, graph.findVertex(0));
    std::cout << "Total distance of TSP tour (Nearest Neighbor): " << totalDistance << " units" << std::endl;
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;

}

void Application::findEulerianCircuit(int u, std::vector<int>& circuit) {
    for (int v : mst[u]) {
        if (!visited[u][v]) {
            visited[u][v] = visited[v][u] = true;
            findEulerianCircuit(v, circuit);
        }
    }
    circuit.push_back(u);
}

void Application::tspChristofides() {
    auto start = std::chrono::steady_clock::now();

    // Step 1: Compute Minimum Spanning Tree (MST)
    primMST();
    std::cout << "prim" << std::endl;
    // Step 2: Find Minimum Weight Perfect Matching of Odd-Degree Vertices
    std::vector<int> oddVertices;
    for (int i = 0; i < distanceMatrix.size(); ++i) {
        if (mst[i].size() % 2 != 0) { // Check for odd degree vertices
            oddVertices.push_back(i);
        }
    }

    std::vector<int> matching(oddVertices.size(), -1); // Initialize matching array
    std::vector<bool> used(oddVertices.size(), false); // Track used vertices in the matching

    // Greedy matching by pairing adjacent odd vertices
    for (int i = 0; i < oddVertices.size(); ++i) {
        if (!used[i]) {
            int closestIndex = -1;
            double minDistance = std::numeric_limits<double>::max();

            // Find the closest unmatched vertex
            for (int j = i + 1; j < oddVertices.size(); ++j) {
                if (!used[j]) {
                    double distance = distanceMatrix[oddVertices[i]][oddVertices[j]];
                    if (distance < minDistance) {
                        minDistance = distance;
                        closestIndex = j;
                    }
                }
            }
            // Pair the vertices
            matching[i] = oddVertices[closestIndex];
            matching[closestIndex] = oddVertices[i];
            used[i] = true;
            used[closestIndex] = true;
        }
    }
    std::cout << "matching" << std::endl;
    // Step 3: Incorporate the matching edges into the MST
    for (int i = 0; i < oddVertices.size(); ++i) {
        if (matching[i] != -1) {
            int u = oddVertices[i];
            int v = matching[i];
            // Add the matching edges to the MST
            mst[u].push_back(v);
            mst[v].push_back(u);
        }
    }
    std::cout << "eulerian" << std::endl;
    // Step 4: Construct an Eulerian Circuit from the MST
    std::vector<int> eulerianCircuit;
    findEulerianCircuit(0, eulerianCircuit);

    // Step 5: Convert the Eulerian Circuit into a Hamiltonian Circuit (TSP Tour)
    std::vector<bool> visited(distanceMatrix.size(), false);
    tspTour.clear();
    for (int vertex : eulerianCircuit) {
        if (!visited[vertex]) {
            tspTour.push_back(vertex);
            visited[vertex] = true;
        }
    }
    tspTour.push_back(eulerianCircuit.front()); // Complete the Hamiltonian circuit
    std::cout << "calculate" << std::endl;
    // Step 6: Calculate the total cost of the TSP Tour
    float totalCost = 0.0f;
    for (size_t i = 0; i < tspTour.size() - 1; ++i) {
        int current = tspTour[i];
        int next = tspTour[i + 1];
        if(distanceMatrix[current][next] > INT_MAX) {
            totalCost += haversineDistance(graph.findVertex(current),graph.findVertex(next));
        }
        else totalCost += distanceMatrix[current][next];
    }


    std::cout << "Total Cost: " << totalCost << std::endl;

    // Measure execution time
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
}

void Application::tspRealWorld(int source) {
    // Check if the graph is Hamiltonian
    std::vector<int> path;
    path.push_back(source);
    int count = 1;
    if (hamiltonianUtil(source, path, visited[source], count)) {
        std::cout << "Graph is Hamiltonian." << std::endl;
        std::cout << "Hamiltonian Cycle: ";
        for (int vertex : path) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "Graph is not Hamiltonian." << std::endl;
    }
}

bool Application::hamiltonianUtil(int v, std::vector<int>& path, std::vector<bool>& visited, int& count) {
    // Base case: If all vertices are included in the path
    if (count == distanceMatrix.size()) {
        // Check if there's an edge from the last vertex in path to the starting vertex
        int startVertex = path.front();
        if (distanceMatrix[v][startVertex] > 0) {
            path.push_back(startVertex);
            return true; // Found a Hamiltonian cycle
        }
        return false;
    }

    // Try different vertices as the next candidate in the Hamiltonian path
    for (int i = 0; i < distanceMatrix.size(); ++i) {
        if (distanceMatrix[v][i] > 0 && !visited[i]) {
            visited[i] = true;
            path.push_back(i);
            count++;

            // Recursively check if this path leads to a Hamiltonian cycle
            if (hamiltonianUtil(i, path, visited, count)) {
                return true;
            }

            // Backtrack if the current vertex doesn't lead to a Hamiltonian cycle
            visited[i] = false;
            path.pop_back();
            count--;
        }
    }

    return false;
}

bool Application::isHamiltonian(const std::vector<int>& path) {
    // A Hamiltonian cycle should contain all vertices exactly once and end at the starting vertex
    if (path.size() != distanceMatrix.size() + 1) {
        return false;
    }

    int startVertex = path.front();
    int endVertex = path.back();

    // Check if the last vertex has an edge to the start vertex to complete the cycle
    return distanceMatrix[endVertex][startVertex] > 0;
}
