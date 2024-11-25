#pragma once
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>

namespace NGnuplot {

namespace {

class TGnuplotPipe {
public:
    explicit TGnuplotPipe(const std::string& command = "gnuplot -persist") 
        : pipe(popen(command.c_str(), "w")) {
        if (!pipe) {
            throw std::runtime_error("Failed to open Gnuplot process.");
        }
    }

    TGnuplotPipe(const TGnuplotPipe&) = delete;
    TGnuplotPipe& operator=(const TGnuplotPipe&) = delete;

    TGnuplotPipe(TGnuplotPipe&& other) noexcept : pipe(other.pipe) {
        other.pipe = nullptr;
    }
    TGnuplotPipe& operator=(TGnuplotPipe&& other) noexcept {
        if (this != &other) {
            if (pipe) {
                pclose(pipe);
            }
            pipe = other.pipe;
            other.pipe = nullptr;
        }
        return *this;
    }

    ~TGnuplotPipe() {
        if (pipe) {
            pclose(pipe);
        }
    }

    void write(const std::string& commands) {
        if (!pipe) {
            throw std::runtime_error("Pipe is closed or invalid.");
        }
        if (fwrite(commands.c_str(), sizeof(char), commands.size(), pipe) != commands.size()) {
            throw std::runtime_error("Failed to write to Gnuplot pipe.");
        }
    }

private:
    FILE* pipe;
};


inline std::string generateGnuplotCommands(const std::string& filename, long double left, long double right) {
    std::ostringstream commands;
    commands << "set term png size 800,600\n";
    commands << "set output '" << filename << ".png'\n";
    commands << "set title 'Displacements and Errors'\n";
    commands << "set xlabel 'Node Position'\n";
    commands << "set ylabel 'Values'\n";
    commands << "set xrange [" << left << ":" << right << "]\n";
    commands << "plot '" << filename << "' using 1:2 with lines title 'Real', "
             << "'" << filename << "' using 1:3 with lines title 'Computed', "
             << "'" << filename << "' using 1:4 with lines title 'Error'\n";
    return commands.str();
}

} // namespace

inline void plot(const std::string& filename, long double left, long double right) {
    std::ifstream dataFile(filename);
    if (!dataFile.is_open()) {
        throw std::runtime_error("Error: Data file '" + filename + "' not found");
    }

    try {
        const std::string commands = generateGnuplotCommands(filename, left, right);

        TGnuplotPipe gnuplot;
        gnuplot.write(commands);
        #ifdef PRINT
        std::cout << "Plot successfully created: " << filename << ".png\n";
        #endif // PRINT
    } catch (const std::exception& error) {
        throw std::runtime_error(std::string("Failed to generate plot: ") + error.what());
    }
}

} // namespace NGnuplot