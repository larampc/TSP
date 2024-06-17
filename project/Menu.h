#ifndef WATERSUPPLYMANAGER_MENU_H
#define WATERSUPPLYMANAGER_MENU_H

/**
 * \mainpage Welcome to the Routing for Ocean Shipping and Urban Deliveries Algorithms Program
 *
 * \section description_sec Project description
 *
 * This program was made for DA course unit of Bachelor in Informatics and Computing Engineering at FEUP.
 * This program has algorithms options for finding optimal or approximated routes for vehicles in generic shipping and delivery scenarios, from urban deliveries to ocean shipping.
 *
 * \section utility_sec What can this project do?
 *
 * Choose the scenario for which you would like to find a route.

 * Run different algorithms such as:
 * - Backtracking;
 * - Triangular Approximation Heuristic;
 * - Nearest Insertion;
 * - Nearest Neighbour;
 * - Other Heuristic;
 * - Two opt.

 * Change settings such as:
 * - Enabling colour mode;
 * - Load a different scenario.
 */


#include "TSP.h"

/**
 * \class Menu
 * \brief This class handles and runs the different commands.
 *
 * This class stores the network and different functions to analyze and get information from the data.
 * It is also responsible for handling inputs and outputs.
 */
class Menu {
private:
    TSP tsp;
public:
    /**
     * \brief Outputs the dataset menu and handles the respective inputs.
     */
    void init();
    /**
     * \brief Outputs the dataset menu for Toy Graphs and handles the respective inputs.
     */
    void toyGraphs();
    /**
     * \brief Outputs the dataset menu for Fully Connected Graphs and handles the respective inputs.
     */
    void fullyConnected();
    /**
     * \brief Outputs the dataset menu for Real World Graphs and handles the respective inputs.
     */
    void realWorld();
    /**
     * \brief Outputs the main menu and handles the respective inputs.
     */
    void run();
    /**
     * \brief Reads an option from the given minimum number to the given maximum from user input.
     *
     * @param low The minimum number.
     * @param n The given number.
     * @return The option from user input.
     */
    static int readOption(int low, int n);
    /**
     * \brief Outputs the settings menu and handles the respective inputs.
     */
    void settings();
    /**
     * \brief Pauses the output until user presses ENTER.
     */
    static void pressEnterToContinue();
    /**
     * \brief Prompts the user to use the two opt algorithm.
     *
     * @return True to use the two opt algorithm, false otherwise.
     */
    bool getTwoOpt();
    /**
     * \brief Prints the given path.
     *
     * @param path The path to print.
     */
    void printPath(const std::vector<unsigned int> &path);
    /**
     * \brief Prompts the user to print the path.
     *
     * @return True to print the path, false otherwise.
     */
    bool seePath();
};


#endif //WATERSUPPLYMANAGER_MENU_H
