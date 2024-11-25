#pragma once
#include <iostream>
#include <optional>

namespace NOptions {

enum class EElementType : int {
    Linear = 1, 
    Cubic = 3 
};

enum class ERestrictionGrade { 
    First, 
    Second, 
    Third 
};

struct TRestriction {
    ERestrictionGrade grade;
    long double position;
    long double value;

    bool operator<(const TRestriction& other) const {
        return position < other.position;
    }
};

namespace {

struct TOptions {
    EElementType type = EElementType::Linear;
    int elementsAmount = 20;
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

inline std::optional<TOptions> getOptions(int argc, char* argv[]) {
    TOptions options;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage();
            options.help = true;
            return options;
        } else if (arg == "-l" || arg == "--linear") {
            options.type = EElementType::Linear;
        } else if (arg == "-c" || arg == "--cubic") {
            options.type = EElementType::Cubic;
        } else if (arg == "-s" || arg == "--size") {
            if (i + 1 < argc) {
                options.elementsAmount = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing value for option -s.\n";
                return std::nullopt;
            }
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return std::nullopt;
        }
    }

    return options;
}

} // namespace NOptions