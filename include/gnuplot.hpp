#pragma once
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>

namespace NGnuplot {

namespace {

class GnuplotPipe {
public:
    explicit GnuplotPipe(const std::string& command = "gnuplot -persist") {
        pipe = popen(command.c_str(), "w");
        if (!pipe) {
            throw std::runtime_error("Failed to open Gnuplot process.");
        }
    }

    GnuplotPipe(const GnuplotPipe&) = delete;
    GnuplotPipe& operator=(const GnuplotPipe&) = delete;

    GnuplotPipe(GnuplotPipe&& other) noexcept : pipe(other.pipe) {
        other.pipe = nullptr;
    }
    GnuplotPipe& operator=(GnuplotPipe&& other) noexcept {
        if (this != &other) {
            if (pipe) {
                pclose(pipe);
            }
            pipe = other.pipe;
            other.pipe = nullptr;
        }
        return *this;
    }

    ~GnuplotPipe() {
        if (pipe) {
            pclose(pipe);
        }
    }

    void write(const std::string& commands) {
        if (pipe) {
            fwrite(commands.c_str(), sizeof(char), commands.size(), pipe);
        } else {
            throw std::runtime_error("Pipe is closed or invalid");
        }
    }

private:
    FILE* pipe;
};

} // namespace

inline void plot(const std::string& filename, double left, double right) {
    std::ifstream dataFile(filename + ".txt");
    if (!dataFile.is_open()) {
        throw std::runtime_error("Error: Data file '" + filename + ".txt' not found");
    }

    try {
        GnuplotPipe gnuplot;

        std::ostringstream commands;
        commands << "set term png size 800,600\n";
        commands << "set output '" << filename << ".png'\n";
        commands << "set title 'Displacements and Errors'\n";
        commands << "set xlabel 'Node Position'\n";
        commands << "set ylabel 'Values'\n";
        commands << "set xrange [" << left << ":" << right << "]\n";
        commands << "plot '" << filename << ".txt' using 1:2 with lines title 'Real', "
                 << "'" << filename << ".txt' using 1:3 with lines title 'Computed', "
                 << "'" << filename << ".txt' using 1:4 with lines title 'Error'\n"
                 ;

        gnuplot.write(commands.str());
        #ifdef PRINT
        std::cout << "Plot successfully created: " << filename << ".png\n";
        #endif // PRINT
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Failed to generate plot: ") + e.what());
    }
}

} // namespace NGnuplot