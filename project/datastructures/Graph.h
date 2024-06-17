#ifndef TSP_DA_GRAPH_H
#define TSP_DA_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <stack>

class Edge;

#define INF std::numeric_limits<double>::max()


/************************* Vertex  **************************/

/**
 * \class Vertex
 * \brief A custom class to represent a Graph's Vertex.
 *
 * This class stores all information and functions of a Graph's Vertex.
 */
class Vertex {
public:
    /**
     * \brief Vertex constructor.
     *
     * @param in The Vertex info.
     */
    explicit Vertex(int in);
    /**
     * \brief Less than operator to compare Vertex dists.
     *
     * @param vertex The Vertex to compare.
     * @return True if this Vertex dist is lower than the dist of the given Vertex.
     */
    bool operator<(Vertex & vertex) const; // required by MutablePriorityQueue

    /**
     * \brief Gets the Vertex info.
     *
     * @return The Vertex info.
     */
    int getInfo() const;
    /**
     * \brief Gets all outgoing Edge from the Vertex.
     *
     * @return All outgoing Edge from the Vertex.
     */
    std::vector<Edge*> getAdj() const;
    /**
     * \brief Gets the Vertex visited state.
     *
     * @return The Vertex visited state.
     */
    bool isVisited() const;
    /**
     * \brief Gets the Vertex processing state.
     *
     * @return The Vertex processing state.
     */
    bool isExplored() const;
    /**
     * \brief Gets the number of incoming Edge to the Vertex.
     *
     * @return The number of incoming Edge to the Vertex.
     */
    unsigned int getIndegree() const;
    /**
     * \brief Gets the Vertex dist.
     *
     * @return The Vertex dist.
     */
    double getDist() const;
    /**
     * \brief Sets the Vertex dist.
     *
     * @param dist The Vertex dist to set.
     */
    void setDist(double dist);
    /**
     * \brief Sets the Vertex longitude.
     *
     * @param longitude The Vertex longitude to set.
     */
    void setLongitude(double longitude);
    /**
     * \brief Sets the Vertex latitude.
     *
     * @param latitude The Vertex latitude to set.
     */
    void setLatitude(double latitude);
    /**
     * \brief Gets the Vertex longitude.
     *
     * @return The Vertex longitude.
     */
    double getLongitude() const;
    /**
     * \brief Gets the Vertex latitude.
     *
     * @return The Vertex latitude.
     */
    double getLatitude() const;
    /**
     * \brief Gets the Vertex path.
     *
     * @return The Vertex path.
     */
    Edge *getPath() const;
    /**
     * \brief Gets all incoming Edge to the Vertex.
     *
     * @return All incoming Edge to the Vertex.
     */
    std::vector<Edge *> getIncoming() const;
    /**
     * \brief Sets the Vertex info.
     *
     * @param info The Vertex info to set.
     */
    void setInfo(int info);
    /**
     * \brief Sets the Vertex visited state.
     *
     * @param visited The Vertex visited state to set.
     */
    void setVisited(bool visited);
    /**
     * \brief Sets the Vertex processing state.
     *
     * @param explored The Vertex processing state to set.
     */
    void setExplored(bool explored);
    /**
     * \brief Sets the Vertex indegree.
     *
     * @param indegree The Vertex indegree to set.
     */
    void setIndegree(unsigned int indegree);
    /**
     * \brief Sets the Vertex path.
     *
     * @param path The Vertex path to set.
     */
    void setPath(Edge *path);
    /**
     * \brief Creates a new Edge from this Vertex to the given Vertex with the given weight.
     *
     * @param dest The destination vertex of the new Edge.
     * @param w The weight of the new Edge.
     * @return The new Edge.
     */
    Edge * addEdge(Vertex *dest, double w);
    /**
     * \brief Deletes all Edge from this Vertex to the Vertex with the given info.
     *
     * @param in The info of the destination Vertex of the Edge to delete.
     * @return True if deleted any Edge, false otherwise.
     */
    bool removeEdge(int in);
    /**
     * \brief Deletes all Edge that are outgoing from this Vertex.
     */
    void removeOutgoingEdges();
    /**
     * \brief Deletes the given Edge from this Vertex.
     *
     * @param edge The Edge to delete.
     */
    void deleteEdge(Edge *edge) const;
    unsigned queueIndex = 0;
protected:
    int info;                // info node
    std::vector<Edge *> adj;  // outgoing edges

    bool visited = false; // used by DFS, BFS, Prim ...
    bool explored = false; // used by isDAG (in addition to the visited attribute)
    unsigned int indegree = 0; // used by topsort
    double dist = 0;
    double latitude = 0, longitude = 0;
    Edge *path = nullptr;

    std::vector<Edge *> incoming; // incoming edges

};

/********************** Edge  ****************************/

/**
 * \class Edge
 * \brief A custom class to represent a Graph's Edge.
 *
 * This class stores all information and functions of a Graph's Edge.
 */
class Edge {
public:
    /**
     * \brief Edge constructor.
     *
     * @param orig The Edge origin Vertex.
     * @param dest The Edge destination Vertex.
     * @param w The Edge weight.
     */
    Edge(Vertex *orig, Vertex *dest, double w);
    /**
     * \brief Gets the Edge destination Vertex.
     *
     * @return The Edge destination Vertex.
     */
    Vertex * getDest() const;
    /**
     * \brief Gets the Edge weight.
     *
     * @return The Edge weight.
     */
    double getWeight() const;
    /**
     * \brief Gets the Edge origin Vertex.
     *
     * @return The Edge origin Vertex.
     */
    Vertex * getOrig() const;
    /**
     * \brief Gets the reverse Edge of this Edge (the Edge that connects the same two Vertex but is in the opposite direction).
     *
     * @return The the reverse Edge of this Edge.
     */
    Edge *getReverse() const;
    /**
     * \brief Sets the reverse Edge of this Edge (the Edge that connects the same two Vertex but is in the opposite direction).
     *
     * @param reverse The the reverse Edge of this Edge to set.
     */
    void setReverse(Edge *reverse);
    /**
     * \brief Sets the Edge weight.
     *
     * @param weight The Edge weight to set.
     */
    void setWeight(double weight);
    /**
     * \brief Gets the Edge active state.
     *
     * @return The Edge active state.
     */
    bool checkActive() const;
    /**
     * \brief Sets the Edge active state.
     *
     * @param active The Edge active state to set.
     */
    void setActive(bool active);
    /**
     * \brief Gets the Edge visited state.
     *
     * @return The Edge visited state.
     */
    bool checkVisited() const;
    /**
     * \brief Sets the Edge visited state to the given state.
     *
     * @param newVisited The state to set the Edge visited state.
     */
    void setVisited(bool newVisited);
protected:
    Vertex * dest; // destination vertex
    double weight; // edge weight, can also be used for capacity

    // used for bidirectional edges
    Vertex *orig;
    Edge *reverse = nullptr;

    bool isActive = true;
    bool visited = false;
};

/********************** Graph  ****************************/

/**
 * \class Graph
 * \brief A custom class to represent a Graph.
 *
 * This class stores all information and functions of a Graph.
 */
class Graph {
public:
    /**
     * \brief Erases all the contents of this graph.
     */
    void reset();
    /**
     * \brief Graph constructor.
     */
    ~Graph();
    /**
     * \brief Auxiliary function to find a Vertex with the given info.
     *
     * @param in The info of the Vertex to find.
     * @return The Vertex with the given info.
     */
    Vertex *findVertex(int in) const;
    /**
     * \brief Auxiliary function to find an Edge with the given origin and destination Vertex.
     *
     * @param source The Edge's origin Vertex.
     * @param dest The Edge's destination Vertex.
     * @return The Edge with the given origin and destination Vertex.
     */
    Edge *findEdge(int source, int dest) const;
    /**
     * \brief Adds a new Vertex with the given info to this Graph.
     *
     * @param in The info of the Vertex to add.
     * @return True if Vertex was added, false otherwise.
     */
    bool addVertex(int in);
    /**
     * \brief Deletes the Vertex with the given info from this Graph.
     *
     * @param in The info of the Vertex to delete.
     * @return True if Vertex was deleted, false otherwise.
     */
    bool removeVertex(int in);
    /**
     * \brief Adds a new Edge with the given weight from the Vertex with the given origin info to the Vertex with the given destination info to this Graph.
     *
     * @param sourc The info of the origin Vertex of the Edge to add.
     * @param dest The info of the destination Vertex of the Edge to add.
     * @param w The weight of the Edge to add.
     * @return True if Edge was added, false otherwise.
     */
    bool addEdge(int sourc, int dest, double w) const;
    /**
     * \brief Adds two new Edge with the given weight from the Vertex with the given origin info to the Vertex with the given destination info and vice-versa to this Graph.
     *
     * @param sourc The info of the origin Vertex of one Edge and destination Vertex of the other Edge to add.
     * @param dest The info of the destination Vertex of one Edge and origin Vertex of the other Edge to add.
     * @param w The weight of the two Edge to add.
     * @return True if both Edge were added, false otherwise.
     */
    bool addBidirectionalEdge(int sourc, int dest, double w) const;
    /**
     * \brief Gets the number of Vertex in this Graph.
     *
     * @return The number of Vertex in this Graph.
     */
    size_t getNumVertex() const;
    /**
     * \brief Gets the VertexSet of this Graph.
     *
     * @return The VertexSet of this Graph.
     */
    std::vector<Vertex*> getVertexSet() const;
    /**
     * \brief Checks if this Graph is a DAG (Directed Acyclic Graph).
     *
     * @return True if this Graph is a DAG, false otherwise.
     */
    bool isDAG() const;
    /**
     * \brief Auxiliary dfs function to check if this Graph is a DAG (Directed Acyclic Graph) from the given starting Vertex.
     *
     * @param v The starting Vertex.
     * @return True if this Graph is a DAG, false otherwise.
     */
    bool dfsIsDAG(Vertex *v) const;
    /**
     * \brief Sorts all Vertex in this Graph in topological order.
     *
     * @return All Vertex in this Graph in topological order.
     */
    std::vector<int> topSort() const;
    /**
     * \brief Sets all Vertex of this Graph visited state to false.
     */
    void resetVisited();
    /**
     * \brief Sets all Edge of this Graph visited state to false.
     */
    void resetEdgeVisited();
    /**
     * \brief Sets all Edge of this Graph active state to true.
     */
    void resetEdgeActive();
    /**
     * \brief Sets all Vertex of this Graph explored state to false.
     */
    void resetExplored();
    /**
     * \brief Performs the Prim's algorithm on this Graph.
     *
     * @return The MST (minimum spanning tree) calculated by the algorithm.
     */
    std::vector<unsigned int> prim();
    /**
     * \brief Initializes the adjacency matrix with zeros and with the given number of Vertex.
     *
     * @param numVertices The number of Vertex to initialize the adjacency matrix.
     */
    void initializeMatrix(int numVertices);
    /**
     * \brief Deletes the adjacency matrix and releases its memory.
     */
    void releaseMemory();
    /**
     * \brief Adds a bidirectional Edge to the adjacency matrix with the given source, destination and weight.
     *
     * @param src The source Vertex of the Edge to add.
     * @param dest The destination Vertex of the Edge to add.
     * @param weight The weight of the Edge to add.
     */
    void matrixAddBidirectionalEdge(int src, int dest, double weight);
    /**
     * \brief Gets the weight of the Edge that has the given source and destination from the adjacency matrix.
     *
     * @param src The source Vertex of the Edge to get the weight from.
     * @param dest The destination Vertex of the Edge to get the weight from.
     * @return The weight of the Edge that has the given source and destination from the adjacency matrix.
     */
    double matrixGetEdgeWeight(int src, int dest);
protected:
    std::vector<Vertex*> vertexSet;    // vertex set
    double** matrix;
};




#endif //TSP_DA_GRAPH_H
