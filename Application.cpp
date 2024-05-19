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
    graphType = i;

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
            Vertex *destVertex = this->graph.findVertex(stoi(saved[1])); //destination

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

    // distance matrix:
    int n = graph.getVertexSet().size();
    visited.resize(n, std::vector<bool>(n, false));
    distanceMatrix.assign(n, std::vector<float>(n, std::numeric_limits<float>::infinity()));

    for(auto vertex: graph.getVertexSet()) {
        distanceMatrix[vertex->getCode()][vertex->getCode()] = 0;
    }
    for(auto vertex: graph.getVertexSet()) {
        for (auto edge: vertex->getAdj()) {
            distanceMatrix[edge->getOrig()->getCode()][edge->getDest()->getCode()] = edge->getWeight();
            distanceMatrix[edge->getDest()->getCode()][edge->getOrig()->getCode()] = edge->getWeight();// reverse

        }
    }
    graph.setDistanceMatrix(distanceMatrix);


}

void Application::resetGraph() {
    for(auto vertex: graph.getVertexSet()) {
        vertex->setVisited(false);
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
    return distance * 1000;
}



//T2.1 .................................................................................................
void Application::tspBacktracking() {

    auto start = std::chrono::steady_clock::now();

    int n = distanceMatrix.size();
    for(auto vertex: graph.getVertexSet())
        vertex->setVisited(false);
    graph.findVertex(0)->setVisited(true);
    float ans = std::numeric_limits<float>::max();

    std::vector<int> path; // to store the current path
    path.push_back(0); // starting from vertex 0

    tspBacktrackingAux(0, n, 1, 0.0f, ans, path);

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Lowest costing path starting at 0: " << ans << std::endl;
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;

}

void Application::tspBacktrackingAux(int currPos, int n, int count, float cost, float& ans, std::vector<int>& path) {

    if (count == n && distanceMatrix[currPos][0] > 0) {
        if (cost + distanceMatrix[currPos][0] < ans) {
            ans = cost + distanceMatrix[currPos][0];
            path.push_back(0); // add start vertex to complete the tour
            path.pop_back(); // remove the start vertex after printing
        }
        return;
    }

    for (auto edge : graph.findVertex(currPos)->getAdj()) {
        int nextVertex = edge->getDest()->getCode();
        if (!edge->getDest()->isVisited() && distanceMatrix[currPos][nextVertex] > 0) {
            edge->getDest()->setVisited(true);
            path.push_back(nextVertex); // add next vertex to the current path
            tspBacktrackingAux(nextVertex, n, count + 1, cost + distanceMatrix[currPos][nextVertex], ans, path);
            path.pop_back(); // backtrack
            edge->getDest()->setVisited(false); // backtrack
        }
    }
}
//T2.2 .................................................................................................
void Application::tspTriangular() { //Complexity – O(V+E) so it is polynomial on the size of G
    if (graphType == 3) {
        std::cout << "This algorithm cannot be performed on the shipping graph, due to the lack of geographical coordinates." << std::endl;
        return;
    }

    resetGraph();
    auto start = std::chrono::steady_clock::now();

    primMST(); // Construct the MST using Prim's algorithm
    std::vector<bool> visited(distanceMatrix.size(), false);
    preorderTraversal(0, visited); // Perform preorder traversal of MST

    if (graphType == 1 || graphType == 2 || graphType == 3) {
        std::cout << "TSP Tour Sequence: ";
        for (size_t i = 0; i < tspTour.size(); i++) {
            std::cout << tspTour[i] << " -> ";
        }
        std::cout << tspTour[0] << std::endl;
    }

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


void Application::primMST() {

    int n = distanceMatrix.size();
    std::vector<bool> mstSet(n, false);
    std::vector<float> key(n, std::numeric_limits<float>::max());
    std::vector<int> parent(n, -1);

    key[0] = 0;
    std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> pq;

    pq.push({0, 0});

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (mstSet[u]) continue;
        mstSet[u] = true;

        for (int v = 0; v < n; v++) {
            if (distanceMatrix[u][v] > 0 && !mstSet[v] && distanceMatrix[u][v] < key[v]) {
                parent[v] = u;
                key[v] = distanceMatrix[u][v];
                pq.push({key[v], v});
            }
        }
    }

    float totalMSTWeight = 0.0f;
    for (int i = 1; i < n; i++) {
        totalMSTWeight += key[i];
    }

    mst.resize(n);
    for (int v = 1; v < n; v++) {
        mst[parent[v]].push_back(v);
        mst[v].push_back(parent[v]);
    }

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


//T2.3 - NEAREST NEIGHBOR................................................................................................................
void Application::tspNearestNeighbor() {
    if (graphType == 3) {
        std::cout << "This algorithm cannot be performed on the shipping graph, due to the lack of geographical coordinates." << std::endl;
        return;
    }
    resetGraph();
    auto start = std::chrono::steady_clock::now();
    int n = distanceMatrix.size();
    int visitedCount = 0;
    Vertex* next;
    Vertex* previous;
    std::vector<bool> visited(n, false);
    double totalDistance = 0.0;

    // start from the first vertex
    Vertex* current = graph.findVertex(0);
    current->setVisited(true);
    visitedCount++;

    while (visitedCount < n) { // visit all nodes:
        bool pathExists = false;
        double minDistance = std::numeric_limits<double>::max();

        for(auto edge: current->getAdj()) { // check adjcent nodes to look for the best local choice
            if(!edge->getDest()->isVisited() && edge->getWeight() < minDistance) {
                next = edge->getDest();
                minDistance = edge->getWeight();
                pathExists = true; // at least one exists!
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
            totalDistance += edge->getWeight();
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

//T2.3 - CHRISTOFIDES................................................................................................................
void Application::tspChristofides() {
    if (graphType == 3) {
        std::cout << "This algorithm cannot be performed on the shipping graph, due to the lack of geographical coordinates." << std::endl;
        return;
    }
    resetGraph();
    auto start = std::chrono::steady_clock::now();

    primMST();

    std::vector<int> oddVertices;
    for (int i = 0; i < distanceMatrix.size(); ++i) {
        if (mst[i].size() % 2 != 0) { // Check for odd degree vertices
            oddVertices.push_back(i);
        }
    }

    // Minimum weight perfect matching of odd-degree vertices
    std::vector<std::pair<int, int>> matchingEdges = blossomAlgorithm(oddVertices);

    // incorporate the matching edges into the MST
    for (const auto &edge : matchingEdges) {
        int u = edge.first;
        int v = edge.second;
        // Add the matching edges to the MST
        mst[u].push_back(v);
        mst[v].push_back(u);
    }

    // construct an Eulerian circuit from the MST
    std::vector<int> eulerianCircuit;
    findEulerianCircuit(0, eulerianCircuit);

    // convert the Eulerian circuit into a Hamiltonian circuit (TSP Tour)
    std::vector<bool> visited(distanceMatrix.size(), false);
    tspTour.clear();
    for (int vertex : eulerianCircuit) {
        if (!visited[vertex]) {
            tspTour.push_back(vertex);
            visited[vertex] = true;
        }
    }

    tspTour.push_back(eulerianCircuit.front()); // complete the Hamiltonian circuit

    // calculate the total cost of the TSP Tour
    float totalCost = 0.0f;
    for (size_t i = 0; i < tspTour.size() - 1; ++i) {
        int current = tspTour[i];
        int next = tspTour[i + 1];
        if(distanceMatrix[current][next] > INT_MAX)
            totalCost += haversineDistance(graph.findVertex(current),graph.findVertex(next));

        else totalCost += distanceMatrix[current][next];
    }

    std::cout << "Total Cost: " << totalCost << std::endl;
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
}

std::vector<std::pair<int, int>> Application::blossomAlgorithm(const std::vector<int>& oddVertices) {
    int n = oddVertices.size(); // Number of odd vertices
    std::vector<std::pair<int, int>> matchingEdges; // Edges of the perfect matching

    // cost -> Distances between the odd vertices
    std::vector<std::vector<double>> cost(n, std::vector<double>(n, std::numeric_limits<double>::max()));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            cost[i][j] = distanceMatrix[oddVertices[i]][oddVertices[j]];
            cost[j][i] = cost[i][j];
        }
    }

    std::vector<int> match(n, -1); // match partner of each odd vertex
    std::vector<bool> used(n, false); // tracks if vertex is already in the matching

    // loop that iterates over the odd vertices until they are all matched
    // iteratively finds augmenting paths

    for (int u = 0; u < n; ++u) {
        if (match[u] == -1) {
            std::vector<double> minCost(n, std::numeric_limits<double>::max());
            std::vector<int> prev(n, -1); // stores previous vertex in the path
            std::vector<bool> inPath(n, false); // tracks if vertex has already been visited in the path

            // double -> cost of reaching a particular vertex
            // int -> index of that vertex
            // the priority queue prioritizes elements with the lowest cost
            std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;

            pq.push({0, u});
            minCost[u] = 0;

            while (!pq.empty()) {
                int v = pq.top().second;
                pq.pop();

                if (inPath[v]) continue;
                inPath[v] = true;

                for (int w = 0; w < n; ++w) {
                    if (!inPath[w] && cost[v][w] < minCost[w]) {
                        minCost[w] = cost[v][w];
                        prev[w] = v;
                        pq.push({minCost[w], w});
                    }
                }
            }

            int v = -1;
            for (int w = 0; w < n; ++w) {
                if (!used[w] && (v == -1 || minCost[w] < minCost[v])) {
                    v = w;
                }
            }

            for (int w = v; w != -1; w = prev[w]) {
                int prev_w = prev[w];
                if (prev_w != -1) {
                    match[w] = prev_w;
                    match[prev_w] = w;
                    used[w] = used[prev_w] = true;
                }
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        if (match[i] != -1 && i < match[i]) {
            matchingEdges.push_back({oddVertices[i], oddVertices[match[i]]});
        }
    }

    return matchingEdges;
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
//T2.4................................................................................................................
void Application::tspRealWorld(int source) {
    resetGraph();
    auto start = std::chrono::steady_clock::now();

    int n = distanceMatrix.size();
    std::vector<int> tour; // store the tour sequence
    tour.push_back(source); // start from source
    graph.findVertex(source)->setVisited(true); // mark source as visited
    int current = source; // initialize the current node to the source

    while (tour.size() < n) { // untill all nodes are visited
        double minDistance = std::numeric_limits<double>::max();
        int nextNode = -1; // initialize the next node

        // iterate through all nodes to find the nearest unvisited node
        for (int i = 0; i < n; ++i) {
            if (!graph.findVertex(i)->isVisited()) {
                double distance = distanceMatrix[current][i];
                if (distance > 0 && distance < minDistance) {
                    minDistance = distance;
                    nextNode = i;
                }
            }
        }

        if (nextNode == -1) {
            std::cout << "The algorithm could not find a Hamiltonian cycle." << std::endl;
            return;
        }

        // mark next node as visited and add it to the tour
        graph.findVertex(nextNode)->setVisited(true);
        tour.push_back(nextNode);
        current = nextNode;
    }

    double totalDistance = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        totalDistance += distanceMatrix[tour[i]][tour[i + 1]];
    }
    totalDistance += distanceMatrix[tour.back()][source];

    std::cout << "TSP Tour Sequence: ";
    for (int node : tour) {
        std::cout << node << " -> ";
    }
    std::cout << source << std::endl;

    std::cout << "Total Distance: " << totalDistance << std::endl;

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;
}
