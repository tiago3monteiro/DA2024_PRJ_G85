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
    std::array<std::string, 7> filenames = {
            "../graphs/tourism.csv",
            "../graphs/stadiums.csv",
            "../graphs/shipping.csv",
            "../graphs/graph1/",
            "../graphs/graph2/",
            "../graphs/graph3/",
            "../graphs/graph4-25/"
    };
    std::cout << i;
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

    /*
    std::cout << "Distance Matrix size:" << graph.getVertexSet().size() << std::endl;
    for(auto i = 0; i < n; i++) {
        for(auto j = 0; j<n;j++) {
            std::cout << distanceMatrix[i][j]<< " ";
        }
        std::cout <<std::endl;
    }
    std::cout << "parsed";*/
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
//execution time 96milissegundos
//grafo 2 -> 529203
//execution time -> 585803 96milissegundos
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
        totalCost += distanceMatrix[current][next];
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
    int n = distanceMatrix.size(); //number of nodes

    std::vector<bool> mstSet(n, false); //visited vertex
    std::vector<float> key(n, std::numeric_limits<float>::max());//lowest edge weight connecting v to a node in the Tree
    std::vector<int> parent(n, -1); //predecessor of v in the Tree

    key[0] = 0; // Start with vertex 0 as the root

    for (int count = 0; count < n - 1; count++) {
        int u = minKey(key, mstSet); //get the closest unprocessed node
        mstSet[u] = true; //visit it

        for (int v = 0; v < n; v++) { // Check if node is not in MST, that exists a direct link between nodes and
            // distance was reduced

            if(distanceMatrix[u][v] >= INT_MAX) {//manually add distance
                Vertex* source = graph.findVertex(u);
                Vertex* destiny = graph.findVertex(v);
                auto distance = haversineDistance(  source, destiny);
                distanceMatrix[u][v] = distance;

            }
            if (distanceMatrix[u][v] > 0 && !mstSet[v] && distanceMatrix[u][v] < key[v]) {
                parent[v] = u; //update parent
                key[v] = distanceMatrix[u][v]; //and change value of lowest edge weight connected to v
            }
        }
    }

    // Construct MST adjacency list for preorder traversal
    mst.resize(n);
    for (int v = 1; v < n; v++) {
        mst[parent[v]].push_back(v);
        mst[v].push_back(parent[v]);
    }
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