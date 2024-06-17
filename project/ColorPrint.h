#ifndef TSP_COLORPRINT_H
#define TSP_COLORPRINT_H

#include <string>

/**
 * \class ColorPrint
 * \brief A custom class to handle prints with colors.
 *
 * This class can print to the console with or without colors.
 */
class ColorPrint {
public:
    static bool colorMode;
    /**
     * \brief Prints the given string with the given color.
     *
     * @param color The color of the letters, if applicable.
     * @param line The string to print.
     */
    ColorPrint(const std::string& color, const std::string& line);
    /**
     * \brief Swaps the printing mode between colored and non colored.
     */
    void static swapColorMode();
};



#endif //TSP_COLORPRINT_H
