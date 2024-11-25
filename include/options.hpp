#pragma once
#include <iostream>
#include <optional>

namespace NOptions {

enum class ElementType : int { 
    Unknown = 0, 
    Linear = 1, 
    Cubic = 3 
};

enum class RestrictGrade { 
    First, 
    Second, 
    Third 
};

struct Restriction {
    RestrictGrade grade;
    long double pos;
    long double val;
};

namespace {

struct Options {
    ElementType type = ElementType::Unknown;
    int elemAmount = 20;
    bool help = false;
};

inline void printBold(const std::string& text) {
    std::cout << "\x1B[1m" << text << "\x1B[0m";
}

inline void printUsage() {
    std::cout << "Usage:\n"
              << "./main [-l | -c] [-s <SIZE>]\n\n"
              << "Options:\n";

    printBold("\t-h, --help");
    std::cout << "\n\t\tdisplay help.\n";

    printBold("\t-s, --size");
    std::cout << "\n\t\tset amount of elements (default is 20).\n";

    printBold("\t-l, --linear");
    std::cout << "\n\t\tUse linear equations for elements (used by default).\n";

    printBold("\t-c, --cubic");
    std::cout << "\n\t\tUse cubic equations for elements.\n";
}

} // namespace

inline std::optional<Options> getOptions(int argc, char* argv[]) {
    Options options;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage();
            options.help = true;
            return options;
        } else if (arg == "-l" || arg == "--linear") {
            options.type = ElementType::Linear;
        } else if (arg == "-c" || arg == "--cubic") {
            options.type = ElementType::Cubic;
        } else if (arg == "-s" || arg == "--size") {
            if (i + 1 < argc) {
                options.elemAmount = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing value for option -s.\n";
                return std::nullopt;
            }
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return std::nullopt;
        }
    }

    if (options.type == ElementType::Unknown) {
        #ifdef PRINT
        std::cout << "Element type not specified. Using default: Linear.\n";
        #endif // PRINT
        options.type = ElementType::Linear;
    }

    return options;
}

} // namespace NOptions