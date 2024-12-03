#pragma once
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>

namespace NOptions {

enum class EElementType : int {
    LINEAR, 
    CUBIC
};

inline std::unordered_map<const EElementType, const std::string> EElementTypeToString {
    { EElementType::LINEAR, "LINEAR" },
    { EElementType::CUBIC, "CUBIC" }
};

inline std::unordered_map<const EElementType, const int> EElementTypeToInt {
    { EElementType::LINEAR, 1 },
    { EElementType::CUBIC, 3 }
};

enum class ERestrictionGrade : int { 
    FIRST, 
    SECOND, 
    THIRD 
};

struct TRestriction {
    ERestrictionGrade Grade;
    long double Position;
    long double Value;

    bool operator<(const TRestriction& other) const {
        return Position < other.Position;
    }
};

namespace {

struct TOptions {
    EElementType Type = EElementType::LINEAR;
    int ElementsAmount = 20;
    bool Help = false;
};

inline void PrintBold(const std::string& text) {
    std::cout << "\x1B[1m" << text << "\x1B[0m";
}

inline void PrintUsage() {
    std::cout << "Usage:\n"
              << "./main [-l | -c] [-s <SIZE>]\n\n"
              << "Options:\n";

    PrintBold("\t-h, --help");
    std::cout << "\n\t\tdisplay help.\n";

    PrintBold("\t-s, --size");
    std::cout << "\n\t\tset amount of elements (default is 20).\n";

    PrintBold("\t-l, --linear");
    std::cout << "\n\t\tUse linear equations for elements (used by default).\n";

    PrintBold("\t-c, --cubic");
    std::cout << "\n\t\tUse cubic equations for elements.\n";
}

} // namespace

inline std::optional<TOptions> GetOptions(int argc, char* argv[]) {
    TOptions options;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            PrintUsage();
            options.Help = true;
            return options;
        } else if (arg == "-l" || arg == "--linear") {
            options.Type = EElementType::LINEAR;
        } else if (arg == "-c" || arg == "--cubic") {
            options.Type = EElementType::CUBIC;
        } else if (arg == "-s" || arg == "--size") {
            if (i + 1 < argc) {
                options.ElementsAmount = std::stoi(argv[++i]);
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