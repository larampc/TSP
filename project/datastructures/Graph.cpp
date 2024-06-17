#include "Graph.h"
#include "MutablePriorityQueue.h"

#include <utility>
#include <stack>

/************************* Vertex  **************************/


Vertex::Vertex(int in): info(in) {}
/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */

Edge * Vertex::addEdge(Vertex *d, double w) {
    auto newEdge = new Edge(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}

/*
 * Auxiliary function to remove an outgoing edge (with a given destination (d))
 * from a vertex (this).
 * Returns true if successful, and false if such edge does not exist.
 */

bool Vertex::removeEdge(int in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge *edge = *it;
        Vertex *dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true; // allows for multiple edges to connect the same pair of vertices (multigraph)
        }
        else {
            it++;
        }
    }
    return removedEdge;
}

/*
 * Auxiliary function to remove an outgoing edge of a vertex.
 */

void Vertex::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge *edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}


bool Vertex::operator<(Vertex & vertex) const {
    return this->dist < vertex.dist;
}

int Vertex::getInfo() const {
    return this->info;
}

std::vector<Edge*> Vertex::getAdj() const {
    return this->adj;
}

bool Vertex::isVisited() const {
    return this->visited;
}

bool Vertex::isExplored() const {
    return this->explored;
}

unsigned int Vertex::getIndegree() const {
    return this->indegree;
}

double Vertex::getDist() const {
    return this->dist;
}

Edge *Vertex::getPath() const {
    return this->path;
}

std::vector<Edge *> Vertex::getIncoming() const {
    return this->incoming;
}

void Vertex::setInfo(int in) {
    this->info = in;
}

void Vertex::setVisited(bool visited) {
    this->visited = visited;
}

void Vertex::setExplored(bool explored) {
    this->explored = explored;
}

void Vertex::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

void Vertex::setDist(double dist) {
    this->dist = dist;
}

void Vertex::setLatitude(double latitude) {
    this->latitude = latitude;
}

void Vertex::setLongitude(double longitude) {
    this->longitude = longitude;
}

void Vertex::setPath(Edge *path) {
    this->path = path;
}

void Vertex::deleteEdge(Edge *edge) const {
    Vertex *dest = edge->getDest();
    // Remove the corresponding edge from the incoming list
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        }
        else {
            it++;
        }
    }
    delete edge;
}

double Vertex::getLongitude() const {
    return longitude;
}

double Vertex::getLatitude() const {
    return latitude;
}

/********************** Edge  ****************************/


Edge::Edge(Vertex *orig, Vertex *dest, double w): orig(orig), dest(dest), weight(w) {}

Vertex * Edge::getDest() const {
    return this->dest;
}

double Edge::getWeight() const {
    return this->weight;
}

Vertex * Edge::getOrig() const {
    return this->orig;
}

Edge *Edge::getReverse() const {
    return this->reverse;
}

void Edge::setReverse(Edge *reverse) {
    this->reverse = reverse;
}

void Edge::setWeight(double weight) {
    this->weight = weight;
}


bool Edge::checkActive() const {
    return isActive;
}

void Edge::setActive(bool active) {
    this->isActive = active;
}

bool Edge::checkVisited() const {
    return visited;
}

void Edge::setVisited(bool newVisited) {
    visited = newVisited;
}

/********************** Graph  ****************************/


size_t Graph::getNumVertex() const {
    return vertexSet.size();
}

std::vector<Vertex*> Graph::getVertexSet() const {
    return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
Vertex* Graph::findVertex(int in) const {
    return (in < vertexSet.size()) ? vertexSet[in] : nullptr;
}


/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
bool Graph::addVertex(int in) {
    if(in < vertexSet.size()) return false;
    vertexSet.push_back(new Vertex(in));
    return true;
}

/*
 *  Removes a vertex with a given content (in) from a graph (this), and
 *  all outgoing and incoming edges.
 *  Returns true if successful, and false if such vertex does not exist.
 */

bool Graph::removeVertex(int in) {
    if (in >= vertexSet.size()) return false;
    Vertex* v = vertexSet[in];
    v->removeOutgoingEdges();
    for (auto u: vertexSet) {
        u->removeEdge(v->getInfo());
    }
    vertexSet.erase(vertexSet.begin() + in);
    delete v;
    return true;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */

bool Graph::addEdge(int sourc, int dest, double w) const {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    v1->addEdge(v2, w);
    return true;
}

bool Graph::addBidirectionalEdge(int sourc, int dest, double w) const {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

/****************** isDAG  ********************/
/*
 * Performs a depth-first search in a graph (this), to determine if the graph
 * is acyclic (acyclic directed graph or DAG).
 * During the search, a cycle is found if an edge connects to a vertex
 * that is being processed in the stack of recursive calls (see theoretical classes).
 * Returns true if the graph is acyclic, and false otherwise.
 */

bool Graph::isDAG() const {
    for (auto v : vertexSet) {
        v->setVisited(false);
        v->setExplored(false);
    }
    return std::none_of(vertexSet.begin(), vertexSet.end(), [this] (Vertex* v){
        return (v->isVisited() && !dfsIsDAG(v));
    });
}

bool Graph::dfsIsDAG(Vertex *v) const {
    v->setVisited(true);
    v->setExplored(true);
    for (auto e : v->getAdj()) {
        if(!e->checkActive()) continue;
        auto w = e->getDest();
        if (w->isExplored()) return false;
        if (! w->isVisited()) {
            if (! dfsIsDAG(w)) return false;
        }
    }
    v->setExplored(false);
    return true;
}

void dfsVisit(Vertex* v, std::stack<int>& aux){
    v->setVisited(true);
    v->setExplored(true);
    for(Edge* adj : v->getAdj()){
        if (adj->checkActive()) {
            if(!adj->getDest()->isVisited()) dfsVisit(adj->getDest(), aux);
        }
    }
    v->setExplored(false);
    aux.push(v->getInfo());
}

std::vector<int> Graph::topSort() const {
    std::vector<int> res;
    std::stack<int> aux;
    for(auto v : vertexSet){
        v->setVisited(false);
        v->setExplored(false);
    }
    for(const auto& v : vertexSet){
        if(!v->isVisited()){
            dfsVisit(v, aux);
        }
    }
    while (!aux.empty()) {
        res.push_back(aux.top());
        aux.pop();
    }
    return res;
}


Edge* Graph::findEdge(int source, int dest) const {
    auto v = findVertex(source);
    for(auto adj: v->getAdj()){
        if(adj->getDest()->getInfo() == dest) return adj;
    }
    return nullptr;
}

void Graph::reset() {
    for (const auto& v: vertexSet) {
        removeVertex(v->getInfo());
    }
    vertexSet.clear();
    releaseMemory();
}

Graph::~Graph() {
    for (const auto& v: vertexSet) {
        removeVertex(v->getInfo());
    }
    vertexSet.clear();
}

std::vector<unsigned int> Graph::prim() {
    MutablePriorityQueue q;
    std::vector<unsigned int> order;
    order.reserve(getNumVertex());
    for(auto v : vertexSet){
        v->setVisited(false);
        v->setDist(INF);
        q.insert(v);
    }
    {
        Vertex *start = findVertex(0);
        start->setDist(0);
        start->setPath(nullptr);
        start->setVisited(true);
        q.decreaseKey(start);
    }
    while(!q.empty()){
        Vertex* u = q.extractMin();
        for(Edge* e : u->getAdj()){
            Vertex* v = e->getDest();
            if(!v->isVisited() && e->getWeight() < v->getDist()){
                v->setPath(e);
                v->setDist(e->getWeight());
                q.decreaseKey(v);
            }
        }
        u->setVisited(true);
        order.push_back(u->getInfo());
    }
    return order;
}

void Graph::resetVisited() {
    for (const auto& v: vertexSet) {
        v->setVisited(false);
    }
}

void Graph::resetEdgeVisited() {
    for (const auto& v: vertexSet) {
        for (const auto &e: v->getAdj()) {
            e->setVisited(false);
        }
    }
}

void Graph::resetEdgeActive() {
    for (const auto& v: vertexSet) {
        for (const auto &e: v->getAdj()) {
            e->setActive(true);
        }
    }
}

void Graph::resetExplored() {
    for (const auto& v: vertexSet) {
        v->setExplored(false);
    }
}

void Graph::initializeMatrix(int numVertices) {
    matrix = new double*[numVertices];
    for (int i = 0; i < numVertices; ++i) {
        matrix[i] = new double[numVertices];
        // Initialize all elements to 0 (assuming no edges initially)
        for (int j = 0; j < numVertices; ++j) {
            matrix[i][j] = -1;
        }
    }
}

void Graph::releaseMemory() {
    for (int i = 0; i < vertexSet.size(); ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void Graph::matrixAddBidirectionalEdge(int src, int dest, double weight) {
    matrix[src][dest] = weight;
    matrix[dest][src] = weight;
}

double Graph::matrixGetEdgeWeight(int src, int dest) {
    return matrix[src][dest];
}

