//
// Created by tiagomonteiro on 5/7/24.
//

#ifndef DA2024_PRJ2_G85__GRAPH_H
#define DA2024_PRJ2_G85__GRAPH_H


#include <iostream>
#include <vector>
#include <queue>

class Edge;

class Vertex {
public:

    Vertex(int c);
    int getCode() const;
    float getLat() const;
    float getLon() const;
    std::vector<Edge *> getAdj();
    bool isVisited() const;
    bool isProcessing() const;
    unsigned int getIndegree() const;
    double getDist() const;
    Edge *getPath() const;
    std::vector<Edge *> getIncoming() const;

    void setCode(int code);
    void setLat(float lat);
    void setLon(float lon);
    void setVisited(bool visited);
    void setProcesssing(bool processing);
    void setIndegree(unsigned int indegree);
    void setDist(double dist);
    void setPath(Edge *path);
    Edge * addEdge(Vertex* dest, double capacity);
    bool removeEdge(int code);
    void removeOutgoingEdges();

    bool operator==(const Vertex& other) const {
        return (this->code == other.code);
    }


protected:
    int code;
    float lat;
    float lon;
    std::vector<Edge*> adj;
    bool visited = false;
    bool processing = false;
    unsigned int indegree;
    double dist = 0;
    Edge *path = nullptr;
    std::vector<Edge *> incoming;
    void deleteEdge(Edge *edge);
};

class Edge {
public:
    Edge(Vertex *orig, Vertex *dest, double capacity);

    Vertex *getDest() const;
    double getCapacity() const;
    bool isSelected() const;
    Vertex * getOrig() const;
    Edge *getReverse() const;
    bool isAnalyzed() const;
    void setAnalyzed(bool analyzed);

protected:
    Vertex * dest;
    double capacity;
    bool selected = false;
    Vertex *orig;
    Edge *reverse = nullptr;
    bool analyzed;
};

class Graph {
public:
   // ~Graph();
    Vertex *findVertex(const int &code) const;
    bool addVertex(Vertex *vertex);
    bool removeVertex(const int code);

    bool addEdge(const int &source, const int &dest, double capacity);
    bool removeEdge(const int &source, const int &dest);
    int getNumVertex() const;
    std::vector<Vertex *> getVertexSet() const;

    bool isDAG() const;
    bool dfsIsDAG(Vertex *v) const;
    std::vector<int> topsort() const;
    void setDistanceMatrix(const std::vector<std::vector<float>> &distMatrix);
    float getDistance(int from, int to) const;

protected:
    std::vector<Vertex *> vertexSet;
    std::vector<std::vector<float>> distanceMatrix;


};


#endif //DA2024_PRJ2_G85__GRAPH_H
