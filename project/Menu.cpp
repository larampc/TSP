#include "Menu.h"
#include <iostream>
#include "ColorPrint.h"

using namespace std;

int Menu::readOption(int low, int n){
    string option;
    getline(cin, option);
    int tmp;
    while( option.empty() || std::any_of(option.begin(),option.end(), [](char c){return !isdigit(c); })
                      || option.size() > to_string(n).size() ||  (tmp = stoi(option)) > n || tmp < low){
        if (!option.empty()) ColorPrint("red", "Invalid Input, please try again\n");
        getline(cin, option);
    }
    return tmp;
}

void Menu::run() {
    while (true) {
        ColorPrint("blue", "\n-----------------------------------\n");
        ColorPrint("blue", "Ocean Shipping and Urban Deliveries Optimizer\n");
        ColorPrint("blue", "-----------------------------------\n");
        ColorPrint("blue", "Select option:\n");
        ColorPrint("cyan", "1. ");
        bool tooBig = tsp.getGraph()->getNumVertex() > 400;
        ColorPrint("white", "TSP Backtracking");
        ColorPrint("yellow", string(tooBig? "(Unadvised - graph too big)":"") + "\n");
        ColorPrint("cyan", "2. ");
        ColorPrint("white", "TSP Backtracking with dp (Held-Karp)");
        ColorPrint("yellow", string(tooBig? "(Unadvised - graph too big)":"") + "\n");
        ColorPrint("cyan", "3. ");
        ColorPrint("white", "Triangular Approximation Heuristic\n");
        ColorPrint("cyan", "4. ");
        ColorPrint("white", "Nearest Insertion\n");
        ColorPrint("cyan", "5. ");
        ColorPrint("white", "Nearest Neighbour\n");
        ColorPrint("cyan", "6. ");
        ColorPrint("white", "Is TSP Possible?\n");
        ColorPrint("cyan", "7. ");
        ColorPrint("white", "Settings \n");
        ColorPrint("cyan", "8. ");
        ColorPrint("red", "Quit Manager \n");
        cin.sync();
        switch (readOption(1, 8)) {
            case 1: {
                std::vector<unsigned int> path(tsp.getGraph()->getNumVertex(), 0);
                double result = tsp.tspBT(path);
                if (result == INF) {
                    ColorPrint("red", "No TSP path found.\n");
                    break;
                }
                ColorPrint("white", "Sum: ");
                ColorPrint("cyan", to_string(result) + "\n");
                if(seePath()) {
                    printPath(path);
                    pressEnterToContinue();
                }
                break;
            }
            case 2: {
                std::vector<unsigned int> path;
                double result = tsp.heldKarp(path);
                if (result == INF) {
                    ColorPrint("red", "No TSP path found.\n");
                    break;
                }
                ColorPrint("white", "Sum: ");
                ColorPrint("cyan", to_string(result) + "\n");
                if(seePath()) {
                    printPath(path);
                    pressEnterToContinue();
                }
                break;
            }
            case 3: {
                bool twoOpt = getTwoOpt();
                std::vector<unsigned int> path;
                double result = tsp.triangularApproximation(path);
                ColorPrint("white", "Result without 2-opt: ");
                ColorPrint("cyan", std::to_string(result) + "\n");
                if (twoOpt) {
                    result = tsp.twoOpt(path);
                    ColorPrint("white", "Result with 2-opt: ");
                    ColorPrint("cyan", std::to_string(result) + "\n");
                }
                if(seePath()) {
                    printPath(path);
                    pressEnterToContinue();
                }
                break;
            }
            case 4: {
                bool twoOpt = getTwoOpt();
                std::vector<unsigned int> path;
                double result = tsp.nearestInsertion(path);
                ColorPrint("white", "Result without 2-opt: ");
                ColorPrint("cyan", std::to_string(result) + "\n");
                if (twoOpt) {
                    result = tsp.twoOpt(path);
                    ColorPrint("white", "Result with 2-opt: ");
                    ColorPrint("cyan", std::to_string(result) + "\n");
                }
                if(seePath()) {
                    printPath(path);
                    pressEnterToContinue();
                }

                break;
            }
            case 5: {
                bool twoOpt = getTwoOpt();
                std::vector<unsigned int> path;
                double result = tsp.nearestNeighbour(path);
                ColorPrint("white", "Result without 2-opt: ");
                ColorPrint("cyan", std::to_string(result) + "\n");
                if (twoOpt) {
                    result = tsp.twoOpt(path);
                    ColorPrint("white", "Result with 2-opt: ");
                    ColorPrint("cyan", std::to_string(result) + "\n");
                }
                if(seePath()) {
                    printPath(path);
                    pressEnterToContinue();
                }

                break;
            }
            case 6: {
                ColorPrint("white", "Enter the vertex you wish to start from \n");
                int start = readOption(0, tsp.getGraph()->getNumVertex());
                vector<unsigned int> path;
                bool result = tsp.realworldTSP(start, path);
                if (result) {
                    if(seePath()) {
                        printPath(path);
                        pressEnterToContinue();
                    }
                }
                else {
                    ColorPrint("red", "No path found! \n");
                    pressEnterToContinue();
                }
                break;
            }
            case 7:
                settings();
                break;
            case 8:
                ColorPrint("blue", "Bye Bye :(\n");
                tsp.getGraph()->releaseMemory();
                exit(0);
        }
    }
}

void Menu::printPath(const std::vector<unsigned int>& path){
    int breakLine = 0;
    ColorPrint("blue", "Path taken: \n");
    ColorPrint("white", to_string(path[0]));
    for (int i = 1; i < tsp.getGraph()->getNumVertex(); i++) {
        ColorPrint("pink", " -> ");
        ColorPrint("white",to_string(path[i]));
        if(++breakLine == 14) {
            std::cout << std::endl;
            breakLine = 0;
        }
    }
    ColorPrint("pink", " -> ");
    ColorPrint("white",to_string(path[0]));
    std::cout << std::endl;
}

bool Menu::getTwoOpt() {
    ColorPrint("blue", "Select option:\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("white", "Without 2-opt \n");
    ColorPrint("cyan", "2. ");
    ColorPrint("white", "With 2-opt\n");
    cin.sync();
    switch (readOption(1, 2)) {
        case 1:
            return false;
        case 2:
            return true;
    }
    return false;
}
bool Menu::seePath() {
    ColorPrint("blue", "See path taken?:\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("white", "Yes \n");
    ColorPrint("cyan", "2. ");
    ColorPrint("white", "No\n");
    cin.sync();
    return (readOption(1,2) == 1);
}

void Menu::settings() {
    string option;
    ColorPrint("blue", "Select option:\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("yellow", ColorPrint::colorMode ? "Disable" : "Enable");
    ColorPrint("white", " Color Mode\n");
    ColorPrint("cyan", "2. ");
    ColorPrint("red", "Load a different graph\n");
    ColorPrint("cyan", "3. ");
    ColorPrint("red", "Cancel\n");
    switch (readOption(1, 3)) {
        case 1:
            ColorPrint::swapColorMode();
            ColorPrint("cyan", ColorPrint::colorMode ? "Color mode enabled\n" : "Color mode disabled\n");
            break;
        case 2:
            init();
            break;
        default:
            break;
    }
}

void Menu::pressEnterToContinue() {
    ColorPrint("cyan", "\nPress ");
    ColorPrint("yellow", "ENTER");
    ColorPrint("cyan", " to continue\n");
    cin.sync();
    cin.ignore();
}


void Menu::init() {
    ColorPrint("blue", "Which data set do you wish to use?\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("white", "Toy Graphs\n");
    ColorPrint("cyan", "2. ");
    ColorPrint("white", "Fully Connected Graphs\n");
    ColorPrint("cyan", "3. ");
    ColorPrint("white", "Real world Graphs\n");
    cin.sync();
    switch (readOption(1, 3)) {
        case 1:
            toyGraphs();
            break;
        case 2:
            fullyConnected();
            break;
        case 3:
            realWorld();
            break;
    }
}

void Menu::toyGraphs() {
    ColorPrint("blue", "Which data set do you wish to use?\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("white", "Shipping\n");
    ColorPrint("cyan", "2. ");
    ColorPrint("white", "Stadiums\n");
    ColorPrint("cyan", "3. ");
    ColorPrint("white", "Tourism\n");
    cin.sync();
    switch (readOption(1,3)) {
        case 1:
            tsp.loadGraph("../dataSets/Toy-Graphs/shipping.csv", true);
            break;
        case 2:
            tsp.loadGraph("../dataSets/Toy-Graphs/stadiums.csv", true);
            break;
        case 3:
            tsp.loadGraph("../dataSets/Toy-Graphs/tourism.csv", true);
            break;
    }
}

void Menu::fullyConnected() {
    ColorPrint("blue", "Which data set do you wish to use?\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("white", "25 edges\n");
    ColorPrint("cyan", "2. ");
    ColorPrint("white", "50 edges\n");
    ColorPrint("cyan", "3. ");
    ColorPrint("white", "75 edges\n");
    ColorPrint("cyan", "4. ");
    ColorPrint("white", "100 edges\n");
    ColorPrint("cyan", "5. ");
    ColorPrint("white", "200 edges\n");
    ColorPrint("cyan", "6. ");
    ColorPrint("white", "300 edges\n");
    ColorPrint("cyan", "7. ");
    ColorPrint("white", "400 edges\n");
    ColorPrint("cyan", "8. ");
    ColorPrint("white", "500 edges\n");
    ColorPrint("cyan", "9. ");
    ColorPrint("white", "600 edges\n");
    ColorPrint("cyan", "10. ");
    ColorPrint("white", "700 edges\n");
    ColorPrint("cyan", "11. ");
    ColorPrint("white", "800 edges\n");
    ColorPrint("cyan", "12. ");
    ColorPrint("white", "900 edges\n");
    cin.sync();
    switch (readOption(1,12)) {
        case 1:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_25.csv", false);
            break;
        case 2:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_50.csv", false);
            break;
        case 3:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_75.csv", false);
            break;
        case 4:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_100.csv", false);
            break;
        case 5:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_200.csv", false);
            break;
        case 6:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_300.csv", false);
            break;
        case 7:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_400.csv", false);
            break;
        case 8:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_500.csv", false);
            break;
        case 9:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_600.csv", false);
            break;
        case 10:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_700.csv", false);
            break;
        case 11:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_800.csv", false);
            break;
        case 12:
            tsp.loadGraph("../dataSets/Extra_Fully_Connected_Graphs/edges_900.csv", false);
            break;
    }
}

void Menu::realWorld() {
    ColorPrint("blue", "Which data set do you wish to use?\n");
    ColorPrint("cyan", "1. ");
    ColorPrint("white", "Graph 1\n");
    ColorPrint("cyan", "2. ");
    ColorPrint("white", "Graph 2\n");
    ColorPrint("cyan", "3. ");
    ColorPrint("white", "Graph 3\n");
    cin.sync();
    switch (readOption(1, 3)) {
        case 1:
            tsp.loadGraphTwoFiles("../dataSets/Real-world Graphs/graph1/nodes.csv", "../dataSets/Real-world Graphs/graph1/edges.csv");
            break;
        case 2:
            tsp.loadGraphTwoFiles("../dataSets/Real-world Graphs/graph2/nodes.csv", "../dataSets/Real-world Graphs/graph2/edges.csv");
            break;
        case 3:
            tsp.loadGraphTwoFiles("../dataSets/Real-world Graphs/graph3/nodes.csv", "../dataSets/Real-world Graphs/graph3/edges.csv");
            break;
    }
}

