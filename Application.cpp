//
// Created by tiagomonteiro on 5/7/24.
//

#include "Application.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <chrono>
#include <array>

Application::Application(int i) {
    std::array<std::string, 6> filenames = {
            "../graphs/tourism.csv",
            "../graphs/stadiums.csv",
            "../graphs/shipping.csv",
            "../graphs/graph1/",
            "../graphs/graph2/",
            "../graphs/graph3/"
    };

    //PARSE TOY GRAPHS:
    if(i == 1 || i == 2 || i == 3)
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
        std::cout << "Real World:" <<std::endl;
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
            vertex->setLon(stof(saved[1]));
            vertex->setLat(stof(saved[2]));
            graph.addVertex(vertex);
        }

        //parse edges:
        std::ifstream in1(edges);
        std::getline(in1,line,'\n'); //remove first line
        while(std::getline(in1,line,'\n')) {

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

    std::cout << "Distance Matrix size:" << graph.getVertexSet().size() << std::endl;



}
//2600 -tourism
//86.7 - shipping
//341 - stadiums


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
    tsp(0, n, 1, 0.0f, ans);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Lowest costing path starting at 0: " << ans << std::endl;
    std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;

}

void Application::tsp(int currPos, int n, int count, float cost, float& ans) {

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
            tsp(nextVertex, n, count + 1, cost + distanceMatrix[currPos][nextVertex], ans);
            edge->getDest()->setVisited(false); //backtrack
        }
    }
}

