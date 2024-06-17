#include <fstream>
#include <sstream>
#include "TSP.h"
#include <cstring>
#include <chrono>
#include <cmath>
#include "ColorPrint.h"

inline static double convert_to_radians(double angle) {
    return angle*M_PI/180;
}

double haversine(double lat1, double lon1, double lat2, double lon2) {
    double rad_lat1 = convert_to_radians(lat1);
    double rad_lon1 = convert_to_radians(lon1);
    double rad_lat2 = convert_to_radians(lat2);
    double rad_lon2 = convert_to_radians(lon2);
    double delta_lat = rad_lat2 - rad_lat1;
    double delta_lon = rad_lon2 - rad_lon1;
    double aux = pow(sin(delta_lat/2), 2) + cos(rad_lat1) * cos(rad_lat2) * pow(sin(delta_lon/2),2);
    double c = 2.0 * atan2( sqrt(aux), sqrt(1.0-aux));
    int earthradius = 6371000;
    return earthradius * c;
}

void TSP::loadGraph(const std::string &path, bool ignoreLine) {
    if(graph.getNumVertex() > 0) graph.reset();
    loadVertexAndEdge(path, ignoreLine);
}
void TSP::loadGraphTwoFiles(const std::string &nodesPath, const std::string &edgesPath) {
    loadVertex(nodesPath);
    loadEdges(edgesPath);
}

void TSP::loadVertex(const std::string &path) {
    std::ifstream file(path);
    std::string line;
    getline(file, line);
    while (getline(file, line)) {
        std::istringstream iss(line);
        int id;
        double latitude, longitude;
        char comma;
        iss >> id >> comma >> latitude >> comma >> longitude;
        graph.addVertex(id);
        Vertex* v = graph.findVertex(id);
        v->setLongitude(longitude);
        v->setLatitude(latitude);
    }
    graph.initializeMatrix(graph.getNumVertex());
    file.close();
}
void TSP::loadEdges(const std::string &path) {
    std::ifstream file(path);
    std::string line;
    getline(file, line);
    while (getline(file, line)) {
        std::istringstream iss(line);
        int id1, id2;
        char comma;
        double dist;
        iss >> id1 >> comma >> id2 >> comma >> dist;
        graph.addBidirectionalEdge(id1, id2, dist);
        graph.matrixAddBidirectionalEdge(id1, id2, dist);
    }
    file.close();
}
void TSP::loadVertexAndEdge(const std::string& path, bool ignoreLine){
    std::ifstream file(path);
    std::string line;
    if (ignoreLine) getline(file, line);
    std::vector<std::tuple<int, int, double>> edges;
    while (getline(file, line)) {
        std::istringstream iss(line);
        int id1, id2;
        char comma;
        double dist;
        iss >> id1 >> comma >> id2 >> comma >> dist;
        edges.emplace_back(id1, id2, dist);
        graph.addVertex(id1);
        if (id2 == graph.getNumVertex()) graph.addVertex(id2);
    }
    graph.initializeMatrix(graph.getNumVertex());
    for (auto e: edges) {
        graph.addBidirectionalEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e));
        graph.matrixAddBidirectionalEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e));
    }
    file.close();
}
double auxBT(unsigned int path[], Vertex* node, double cost, unsigned int k, unsigned int n, double minimum) {
    if (cost > minimum) return INF;
    if (k == n ) {
        for (auto e: node->getAdj()) {
            if (e->getDest()->getInfo() == 0) return cost + e->getWeight();
        }
        return INF;
    }
    unsigned int bestPath[n];
    for (auto e: node->getAdj()) {
        if (!e->getDest()->isVisited()) {
            e->getDest()->setVisited(true);
            unsigned int tmp[n];
            memcpy(tmp, path, n*sizeof(unsigned int));
            tmp[k] = e->getDest()->getInfo();
            double value = auxBT(tmp, e->getDest(), cost + e->getWeight(), k+1, n, minimum);
            if (value < minimum) {
                memcpy(bestPath, tmp, n*sizeof (unsigned int));
                minimum = value;
            }
            e->getDest()->setVisited(false);
        }
    }
    memcpy(path, bestPath, n*sizeof(unsigned int));
    return minimum;
}


double TSP::tspBT(std::vector<unsigned int>& path) {
    unsigned int tour[path.size()];

    auto ti = std::chrono::high_resolution_clock::now();
    memset(tour, 0, graph.getNumVertex() * sizeof (unsigned int));
    graph.resetVisited();
    Vertex* start = graph.findVertex(0);
    start->setVisited(true);
    double minimum = INF;
    double best = auxBT(tour, start, 0, 1, graph.getNumVertex(), minimum);
    auto tf = std::chrono::high_resolution_clock::now();

    ColorPrint("white", "Execution time: ");
    ColorPrint("cyan", std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tf-ti).count()) + " ms\n");
    for(int i = 0; i < path.size(); i++) path[i] = tour[i];
    return best;
}

Graph* TSP::getGraph() {
    return &graph;
}

double TSP::triangularApproximation(std::vector<unsigned int>& path){
    auto ti = std::chrono::high_resolution_clock::now();
    path = graph.prim();
    path.push_back(0);
    double sum = 0;
    bool usedHaversine = false;
    for(int i = 0; i < path.size()-1; i++){
        sum += getDistMatrix(path[i], path[i+1], usedHaversine);
    }
    auto tf = std::chrono::high_resolution_clock::now();
    if (usedHaversine) ColorPrint("yellow", "Warning: Used Haversine to compute distance between non existing edges. \n");
    ColorPrint("white", "Execution time: ");
    ColorPrint("cyan", std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tf-ti).count()) + " ms\n");
    return sum;
}

void TSP::do2Opt(std::vector<unsigned int>& tour, int i, int j){
    std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
}

double TSP::dist(unsigned int v1, unsigned int v2){
    Edge* e = graph.findEdge(v1,v2);
    Vertex* u = graph.findVertex(v1);
    Vertex* v = graph.findVertex(v2);
    return e != nullptr ? e->getWeight() : haversine(u->getLatitude(), u->getLongitude(), v->getLatitude(), v->getLongitude());
}

double TSP::twoOpt(std::vector<unsigned int>& tour) {
    double improvement = 0;
    bool improved = true;
    int tries = 10;
    size_t n = graph.getNumVertex();
    while (improved && tries--) {
        improved = false;
        for (int i = 0; i <= n - 2; i++) {
            for (int j = i + 1; j <= n - 1; j++) {
                double lengthDelta = -dist(tour[i],tour[i + 1]) - dist(tour[j],tour[j + 1])
                                     + dist(tour[i],tour[j]) + dist(tour[i + 1],tour[j + 1]);

                if (lengthDelta < 0) {
                    do2Opt(tour, i, j);
                    improvement += lengthDelta;
                    improved = true;
                }
            }
        }
    }
    double sum = 0;
    bool usedHaversine = false;
    for(int i = 0; i < tour.size()-1; i++){
        sum += getDistMatrix(tour[i], tour[i+1], usedHaversine);
    }
    if (usedHaversine) ColorPrint("yellow", "Warning: Used Haversine to compute distance between non existing edges. \n");
    return sum;
}

double TSP::getDistMatrix(int src, int dest, bool& usedHaversine) {
    if (graph.matrixGetEdgeWeight(src, dest) >= 0) return graph.matrixGetEdgeWeight(src, dest);
    else {
        usedHaversine = true;
        return haversine(graph.findVertex(src)->getLatitude(), graph.findVertex(src)->getLongitude(), graph.findVertex(dest)->getLatitude(), graph.findVertex(dest)->getLongitude());
    }
}

double TSP::nearestInsertion(std::vector<unsigned int>& path) {
    auto ti = std::chrono::high_resolution_clock::now();
    std::vector<Edge*> edges;
    for (const auto v: graph.getVertexSet()) {
        v->setIndegree(0);
        v->setVisited(false);
        for (auto e: v->getAdj()) {
            e->setVisited(false);
        }
    }
    path.push_back(0);
    graph.findVertex(0)->setVisited(true);
    double total = 0;
    bool usedHaversine = false;
    while (path.size() < graph.getNumVertex()) {
        double min = INF;
        std::pair<int, int> best = {};
        for (auto node: path) {
            Vertex* v = graph.findVertex(node);
            if (v->getIndegree() >= 2) continue;
            for (int i = 0; i < graph.getNumVertex(); i++) {
                if (i == v->getInfo()) continue;
                bool preUsedHaversine = false;
                double weight = getDistMatrix(v->getInfo(), i, preUsedHaversine);
                if (weight < min && !graph.findVertex(i)->isVisited() && v->getIndegree() < 2) {
                    if (!usedHaversine) usedHaversine = preUsedHaversine;
                    min = weight;
                    best.first = v->getInfo();
                    best.second = i;
                }
            }
        }
        total += min;
        graph.findVertex(best.second)->setIndegree(graph.findVertex(best.second)->getIndegree()+1);
        graph.findVertex(best.first)->setIndegree(graph.findVertex(best.first)->getIndegree()+1);
        graph.findVertex(best.second)->setVisited(true);
        path.push_back(best.second);
    }
    for (auto node: path) {
        Vertex* v = graph.findVertex(node);
        if (v->getIndegree() < 2) {
            for (int i = 0; i < graph.getNumVertex(); i++) {
                if (i == graph.findVertex(*(path.end()-1))->getInfo()) {
                    total += getDistMatrix(v->getInfo(), i, usedHaversine);
                }
            }
            break;
        }
    }
    path.push_back(0);
    auto tf = std::chrono::high_resolution_clock::now();
    if (usedHaversine) ColorPrint("yellow", "Warning: Used Haversine to compute distance between non existing edges. \n");
    ColorPrint("white", "Execution time: ");
    ColorPrint("cyan", std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tf-ti).count()) + " ms\n");
    return total;
}

double TSP::nearestNeighbour(std::vector<unsigned int>& path) {
    auto ti = std::chrono::high_resolution_clock::now();
    graph.resetVisited();
    Vertex* start = graph.findVertex(0);
    start->setVisited(true);
    path.push_back(0);
    int count  = 1;
    double total = 0;
    bool usedHaversine = false;
    while (count < graph.getNumVertex()) {
        double min = INF;
        Vertex* next = nullptr;
        for (int i = 0; i < graph.getNumVertex(); i++) {
            if (i == start->getInfo()) continue;
            bool preUsedHaversine = false;
            double weight = getDistMatrix(start->getInfo(), i, preUsedHaversine);
            if (weight < min && !graph.findVertex(i)->isVisited()) {
                if (!usedHaversine) usedHaversine = preUsedHaversine;
                min = weight;
                next = graph.findVertex(i);
            }
        }
        next->setVisited(true);
        start = next;
        total += min;
        path.push_back(next->getInfo());
        count++;
    }
    total += getDistMatrix(start->getInfo(), 0, usedHaversine);
    path.push_back(0);
    auto tf = std::chrono::high_resolution_clock::now();
    if (usedHaversine) ColorPrint("yellow", "Warning: Used Haversine to compute distance between non existing edges. \n");
    ColorPrint("white", "Execution time: ");
    ColorPrint("cyan", std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tf-ti).count()) + " ms\n");
    return total;
}

bool compareEdges(Edge *e1, Edge *e2)
{
    return (e1->getWeight() < e2->getWeight());
}

std::vector<Vertex *> TSP::realworldStartSet(int start, int end, double& sum) {
    std::vector<Vertex *> res;
    graph.findVertex(start)->setVisited(true);
    std::vector<Edge *> adj = graph.findVertex(start)->getAdj();
    std::sort(adj.begin(), adj.end(), compareEdges);
    for (const auto &e: adj) {
        if (e->getDest()->getInfo() == end && sum > 1) {
            res.push_back(graph.findVertex(start));
            e->setVisited(true);
            e->getReverse()->setVisited(true);
            sum = e->getWeight();
            return res;
        }
        if (!e->getDest()->isVisited()) {
            sum += 1;
            res = realworldStartSet(e->getDest()->getInfo(), end, sum);
            if (!res.empty()) {
                res.push_back(graph.findVertex(start));
                e->setVisited(true);
                e->getReverse()->setVisited(true);
                sum += e->getWeight();
                return res;
            }
            sum -= 1;
        }
    }
    return {};
}

std::vector<Vertex *> TSP::realworldExpandSet(int realStart, int start, std::unordered_set<int> &cycle, std::vector<int> &end, double &sum) {
    std::vector<Vertex *> res;
    graph.findVertex(start)->setVisited(true);
    std::vector<Edge *> adj = graph.findVertex(start)->getAdj();
    std::sort(adj.begin(), adj.end(), compareEdges);
    for (const auto &e: adj) {
        if (realStart != start && std::find(end.begin(), end.end(), e->getDest()->getInfo()) != end.end()) {
            res.push_back(graph.findVertex(start));
            e->setVisited(true);
            e->getReverse()->setVisited(true);
            sum += e->getWeight();
            end = {e->getDest()->getInfo()};
            return res;
        }
        if (!e->getDest()->isVisited() && !e->getDest()->isExplored() && std::find(end.begin(), end.end(), e->getDest()->getInfo()) == end.end() && !cycle.count(e->getDest()->getInfo())) {
            res = realworldExpandSet(realStart, e->getDest()->getInfo(), cycle, end, sum);
            if (!res.empty()) {
                res.push_back(graph.findVertex(start));
                e->setVisited(true);
                e->getReverse()->setVisited(true);
                sum += e->getWeight();
                return res;
            }
        }
    }
    return {};
}

std::vector<Vertex *> TSP::realworldSwitchPath(int realStart, int start, std::string path, std::unordered_set<std::string> &donePaths, std::unordered_set<int> &cycle, double &sum) {
    std::vector<Vertex *> res;
    graph.findVertex(start)->setVisited(true);
    std::vector<Edge *> adj = graph.findVertex(start)->getAdj();
    std::sort(adj.begin(), adj.end(), compareEdges);
    for (const auto &e: adj) {
        if (realStart != start && !e->getDest()->isVisited() && e->getDest()->isExplored() && !donePaths.count(path + std::to_string(e->getDest()->getInfo()))) {
            res.push_back(graph.findVertex(start));
            Vertex* v1 = nullptr;
            Vertex* v2 = nullptr;
            for (const auto &ee: e->getDest()->getAdj()) {
                if (ee->checkVisited()){
                    if (v1 == nullptr) v1 = ee->getDest();
                    else {
                        v2 = ee->getDest();
                        break;
                    }
                }
            }
            Vertex* prev1 = e->getDest();
            Vertex* prev2 = e->getDest();
            while (prev1->getInfo() != realStart && prev2->getInfo() != realStart) {
                for (const auto &ee: v1->getAdj()) {
                    if (ee->checkVisited() && ee->getDest() != prev1){
                        prev1 = v1;
                        v1 = ee->getDest();
                        break;
                    }
                }
                for (const auto &ee: v2->getAdj()) {
                    if (ee->checkVisited() && ee->getDest() != prev2){
                        prev2 = v2;
                        v2 = ee->getDest();
                        break;
                    }
                }
            }
            if (prev1->getInfo() != realStart) {
                v1 = prev2;
                prev1 = v2;
            } else {
                v2 = v1;
                v1 = prev1;
                prev1 = v2;
            }
            std::string path1;
            std::string path2;
            while (v1 != e->getDest()) {
                path1 += std::to_string(v1->getInfo());
                path2 = std::to_string(v1->getInfo()) + path2;
                for (const auto &ee: v1->getAdj()) {
                    if (ee->checkVisited() && ee->getDest() != prev1){
                        prev1 = v1;
                        v1 = ee->getDest();
                        sum -= ee->getWeight();
                        ee->setVisited(false);
                        ee->getReverse()->setVisited(false);
                        if (v1 != e->getDest()) {
                            cycle.erase(v1->getInfo());
                            v1->setExplored(false);
                        }
                        break;
                    }
                }
            }
            path1 += std::to_string(e->getDest()->getInfo());
            path2 = std::to_string(e->getDest()->getInfo()) + path2;
            donePaths.emplace(path1);
            donePaths.emplace(path2);
            donePaths.emplace(path + std::to_string(e->getDest()->getInfo()));

            e->setVisited(true);
            e->getReverse()->setVisited(true);
            sum += e->getWeight();
            return res;
        }
        if (!e->getDest()->isVisited() && !e->getDest()->isExplored()) {
            res = realworldSwitchPath(realStart, e->getDest()->getInfo(), path + std::to_string(e->getDest()->getInfo()), donePaths, cycle, sum);
            if (!res.empty()) {
                res.push_back(graph.findVertex(start));
                e->setVisited(true);
                e->getReverse()->setVisited(true);
                sum += e->getWeight();
                return res;
            }
        }
    }
    return {};
}

bool TSP::realworldTSP(int start, std::vector<unsigned int> &path) {
    auto ti = std::chrono::high_resolution_clock::now();

    double sum = 0;
    int realStart = start;
    graph.resetExplored();
    graph.resetVisited();
    graph.resetEdgeVisited();
    std::unordered_set<int> cycle;
    std::unordered_set<std::string> donePaths;
    cycle.emplace(start);
    std::vector<Vertex *> res = realworldStartSet(start, start, sum);
    if (res.empty()) return false;
    for (const auto &v: res) {
        cycle.emplace(v->getInfo());
    }
    for (const auto &v: res) {
        if (!v->isExplored()) {
            bool explored = true;
            for (Edge *e: v->getAdj()) {
                if (!cycle.count(e->getDest()->getInfo())) {
                    explored = false;
                    break;
                }
            }
            v->setExplored(explored);
        }
    }
    unsigned int total = graph.getNumVertex();
    while (cycle.size() != total) {
        start = -1;
        for (const auto &v: cycle) {
            if (!graph.findVertex(v)->isExplored()) {
                bool valid = false;
                for (const auto &e: graph.findVertex(v)->getAdj()) {
                    if (e->checkVisited()) {
                        if (!e->getDest()->isExplored()) {
                            valid = true;
                            break;
                        }
                    }
                }
                if (valid) {
                    start = v;
                    break;
                } else graph.findVertex(v)->setExplored(true);
            }
        }
        if (start == -1) {
            bool end = true;
            for (const auto &v: cycle) {
                bool valid = false;
                for (const auto &e: graph.findVertex(v)->getAdj()) {
                    if (!e->getDest()->isExplored()) {
                        valid = true;
                        break;
                    }
                }
                if (valid) {
                    graph.resetVisited();
                    res = realworldSwitchPath(v, v, std::to_string(v), donePaths, cycle, sum);
                    if (!res.empty()) {
                        res.pop_back();
                        for (const auto &vv: res) {
                            cycle.emplace(vv->getInfo());
                        }
                        for (const auto &vv: cycle) {
                            bool explored = true;
                            for (Edge *e: graph.findVertex(vv)->getAdj()) {
                                if (!cycle.count(e->getDest()->getInfo())) {
                                    explored = false;
                                    break;
                                }
                            }
                            graph.findVertex(vv)->setExplored(explored);
                        }
                        end = false;
                        break;
                    }
                } else graph.findVertex(v)->setExplored(true);
            }
            if (end) return false;
            else continue;
        }

        std::vector<int> end;
        for (const auto &e: graph.findVertex(start)->getAdj()) {
            if (e->checkVisited() && !e->getDest()->isExplored()) {
                end.push_back(e->getDest()->getInfo());
            }
        }
        graph.resetVisited();
        res = realworldExpandSet(start, start, cycle, end, sum);
        if (res.empty()) {
            graph.findVertex(start)->setExplored(true);
            continue;
        }
        res.pop_back();
        graph.findEdge(start, *end.begin())->setVisited(false);
        graph.findEdge(*end.begin(), start)->setVisited(false);
        sum -= graph.findEdge(start, *end.begin())->getWeight();
        auto endV = *end.begin();
        for (const auto &v: res) {
            cycle.emplace(v->getInfo());
        }
        res.push_back(graph.findVertex(start));
        res.push_back(graph.findVertex(endV));
        for (const auto &v: res) {
            if (!v->isExplored()) {
                bool explored = true;
                for (Edge *e: v->getAdj()) {
                    if (!cycle.count(e->getDest()->getInfo())) {
                        explored = false;
                        break;
                    }
                }
                v->setExplored(explored);
            }
        }
    }
    auto tf = std::chrono::high_resolution_clock::now();
    ColorPrint("white", "Execution time: ");
    ColorPrint("cyan", std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tf-ti).count()) + " ms\n");
    ColorPrint("white", "Sum: ");
    ColorPrint("cyan", std::to_string(sum) + "\n");

    path.push_back(realStart);
    Vertex* v = nullptr;
    Vertex* prev = graph.findVertex(realStart);
    for (const auto &e: prev->getAdj()) {
        if (e->checkVisited()){
            v = e->getDest();
            break;
        }
    }
    while (v->getInfo() != realStart) {
        for (const auto &e: v->getAdj()) {
            if (e->checkVisited() && e->getDest() != prev){
                path.push_back(v->getInfo());
                prev = v;
                v = e->getDest();
                break;
            }
        }
    }
    return true;
}


double TSP::heldKarp(std::vector<unsigned int>& tour)
{
    auto ti = std::chrono::high_resolution_clock::now();

    size_t n = graph.getNumVertex();
    std::vector<std::vector<double>> memo(n, std::vector<double>(1 << n, -1));
    std::vector<std::vector<int>> path = std::vector<std::vector<int>>(n, std::vector<int>(1 << n, -1));
    double result = heldKarp(0, 1, memo, path);
    auto tf = std::chrono::high_resolution_clock::now();


    int mask = 1;
    int k = 0;
    while (true) {
        tour.push_back(k);
        if (mask == (1 << n) - 1) break;
        k = path[k][mask];
        mask |= (1 << k);
    }

    tour.push_back(0);

    ColorPrint("white", "Execution time: ");
    ColorPrint("cyan", std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(tf-ti).count()) + " ms\n");
    return result;
}

double TSP::heldKarp(int curr, unsigned long long mask, std::vector<std::vector<double>> &memo, std::vector<std::vector<int>>& path)
{
    if (memo[curr][mask] != -1) return memo[curr][mask];

    double minimum = INF, tmp;

    if (mask == (1 << graph.getNumVertex()) - 1) {
        return graph.matrixGetEdgeWeight(curr, 0) >= 0 ? graph.matrixGetEdgeWeight(curr, 0) : INF;
    }

    for (int i = 0; i < graph.getNumVertex(); i++) {
        if (graph.matrixGetEdgeWeight(curr, i) < 0) continue;
        if (!(mask & (1 << i))) {
            tmp = graph.matrixGetEdgeWeight(curr, i) + heldKarp(i, mask | (1 << i), memo, path);
            if (tmp < minimum) {
                minimum = tmp;
                path[curr][mask] = i;
            }
        }
    }

    return (memo[curr][mask] = minimum);
}
