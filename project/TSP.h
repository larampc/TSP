#ifndef TSP_DA_TSP_H
#define TSP_DA_TSP_H


#include "Graph.h"

/**
 * \class TSP
 * \brief This class includes several functions that perform algorithms for finding the desired routes.
 *
 * This class stores the graph from the chosen scenario and different functions to obtain a route, solving in a optimal or approximated way the Traveling Salesman Problem.
 */
class TSP {
private:
    Graph graph;
    /**
     * \brief Reads the Vertex file from the given filepath and adds them to the Graph.
     *
     * @param path The filepath of the Vertex file to read from.
     */
    void loadVertex(const std::string& path);
    /**
     * \brief Reads the Edge file from the given filepath and adds them to the Graph.
     *
     * @param path The filepath of the Edge file to read from.
     */
    void loadEdges(const std::string &path);
    /**
     * \brief Reads the Vertex and Edge file from the given filepath and adds them to the Graph.
     *
     * @param path The filepath of the Vertex and Edge file to read from.
     * @param ignoreLine True if the file's first line should be ignored, false otherwise.
     */
    void loadVertexAndEdge(const std::string &path, bool ignoreLine);
public:
    /**
     * \brief Loads the graph which edges are given in the path file.
     *
     * @param path The file path to load the edges from.
     * @param ignoreLine True if the file's first line should be ignored, false otherwise.
     */
    void loadGraph(const std::string& path, bool ignoreLine);
    /**
     * \brief Loads the graph which nodes are given in the nodesPath file and edges in the edgesPath file.
     *
     * @param nodesPath The file path to load the nodes from.
     * @param edgesPath The file path to load the edges from.
     */
    void loadGraphTwoFiles(const std::string &nodesPath, const std::string &edgesPath);
    /**
     * \brief Gets the current graph.
     *
     * @return The current Graph.
     */
    Graph* getGraph();
    /**
     * \brief Performs the optimal Backtracking algorithm to the TSP for a graph starting and ending the node tour on node labelled with the zero-identifier label.
     *
     * @param path The vector where the found path will be stored.
     * @return The final weight of the given path.
     *
     * \par Complexity
     * O(V!) in which V is the number of vertexes.
     */
    double tspBT(std::vector<unsigned int> &path);
    /**
     * \brief Approximates a solution to the TSP using the Nearest Insertion Heuristic.
     *
     * @param path The vector where the found path will be stored.
     * @return The weight of the resulting path.
     *
     * \par Complexity
     * O(V^2) in which V is the number of vertexes.
     */
    double nearestInsertion(std::vector<unsigned int> &path);
    /**
     * \brief Approximates a solution to the TSP using the Nearest Neighbour Heuristic.
     *
     * @param path The vector where the found path will be stored.
     * @return The weight of the resulting path.
     *
     * \par Complexity
     * O(V^2) in which V is the number of vertexes.
     */
    double nearestNeighbour(std::vector<unsigned int> &path);
    /**
     * \brief Auxiliary recursive function to start a cycle from the given start Vertex to a given possible end Vertex.
     *
     * @param start The Vertex to start the cycle from.
     * @param end The possible end Vertex to end the cycle.
     * @param sum The sum of the created cycle.
     * @return The created cycle.
     *
     * \par Complexity
     * O((V log V) * E) in which V is the number of vertexes and E is the number of edges.
     */
    std::vector<Vertex *> realworldStartSet(int start, int end, double& sum);
    /**
     * \brief Auxiliary recursive function to expand the current cycle from the given start Vertex to one of the given end Vertex.
     *
     * @param realStart The Vertex to start expanding from.
     * @param start The current start Vertex, (can change each recursive call).
     * @param cycle The cycle to expand.
     * @param end The possible end Vertex to end the expansion.
     * @param sum The sum of the cycle to expand.
     * @return The expansion made to the cycle.
     *
     * \par Complexity
     * O((V log V) * E) in which V is the number of vertexes and E is the number of edges.
     */
    std::vector<Vertex *> realworldExpandSet(int realStart, int start, std::unordered_set<int> &cycle, std::vector<int> &end, double& sum);
     /**
      * \brief Auxiliary recursive function to modify the current cycle from the given start Vertex.
      *
      * @param realStart The Vertex to start modifying from.
      * @param start The current start Vertex, (can change each recursive call).
      * @param path The current modification to the cycle, (can change each recursive call).
      * @param donePaths Prohibited paths, (path cannot be equal to any of these paths)
      * @param cycle The cycle to modify.
      * @param sum The sum of the cycle to modify.
      * @return The modification made to the cycle.
      *
      * \par Complexity
      * O(E^2 * V) in which V is the number of vertexes and E is the number of edges.
      */
    std::vector<Vertex *> realworldSwitchPath(int realStart, int start, std::string path, std::unordered_set<std::string> &donePaths, std::unordered_set<int> &cycle, double &sum);
    /**
     * \brief Approximates a solution to the TSP starting at the given Vertex returning whether or not it is possible.
     *
     * @param start The Vertex to start the algorithm from.
     * @param path The resulting path, (empty if not possible).
     * @return True if it is possible to compute a path on the Graph false otherwise.
     *
     * \par Complexity
     * O(E^2 * V^3) in which V is the number of vertexes and E is the number of edges.
     */
    bool realworldTSP(int start, std::vector<unsigned int> &path);
    /**
     * \brief Returns the distance between the given source and destination Vertex using the distance matrix. If an edge between them does not exist, computes it via the geographic Haversine distance between them.
     *
     * @param src The source Vertex to get the distance.
     * @param dest The destination Vertex to get the distance.
     * @param usedHaversine True if used Haversine false otherwise.
     * @return The distance between the given source and destination Vertex.
     */
    double getDistMatrix(int src, int dest, bool& usedHaversine);
    /**
     * \brief Approximates a solution to the TSP using the Triangular Approximation Heuristic for a graph starting and ending the node tour on node labelled with the zero-identifier label. Relies on the triangular inequality.
     *
     * @param path The vector where the found path will be stored.
     * @return The sum of the computed path.
     *
     * \par Complexity
     * O(E log V) in which V is the number of vertexes and E the number of edges.
     */
    double triangularApproximation(std::vector<unsigned int> &path);
    /**
     * \brief Performs an optimization of the given tour using the 2-opt algorithm.
     *
     * @param tour The tour to be optimized.
     * @return The sum of the optimised tour.
     *
     * \par Complexity
     * O(V^2) in which V is the number of vertexes.
     */
    double twoOpt(std::vector<unsigned int> &tour);
    /**
     * \brief Returns the distance between the given source and destination Vertex. If an edge between them does not exist, computes it via the geographic Haversine distance between them.
     *
     * @param v1 The source vertex of the edge.
     * @param v2 The destination vertex of the edge.
     * @return The distance between the given vertices.
     */
    double dist(unsigned int v1, unsigned int v2);
    /**
     * \brief Performs a 2-opt swap, swapping the vertices i and j and manipulating the tour as needed.
     *
     * @param tour The tour to be manipulated.
     * @param i The vertex of the first edge to be swapped.
     * @param j The vertex of the second edge to be swapped.
     */
    void do2Opt(std::vector<unsigned int> &tour, int i, int j);
    /**
     * \brief Performs the Held-Karp dynamic programming algorithm to the TSP for a graph starting and ending the node tour on node labelled with the zero-identifier label.
     *
     * @param path The vector where the found path will be stored.
     * @return The final weight of the given path.
     *
     * \par Complexity
     * O(V^2 * 2^V) in which V is the number of vertexes.
     */
    double heldKarp(std::vector<unsigned int> &tour);
    /**
     * \brief Auxiliary Recursive function for performing the Held-Karp dynamic programming algorithm to the TSP for a graph starting and ending the node tour on node labelled with the zero-identifier label.
     *
     * @param curr The current vertex being explored.
     * @param mask The mask storing which vertices have been explored.
     * @param memo The memoization table storing calculated results.
     * @param path The vector with the current path.
     * @return The step weight of the given state.
     */
    double heldKarp(int curr, unsigned long long int mask, std::vector<std::vector<double>> &memo,
                    std::vector<std::vector<int>> &path);
};


#endif //TSP_DA_TSP_H
