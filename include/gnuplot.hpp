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
        : Pipe_(popen(command.c_str(), "w")) {
        if (!Pipe_) {
            throw std::runtime_error("Failed to open Gnuplot process.");
        }
    }

    TGnuplotPipe(const TGnuplotPipe&) = delete;
    TGnuplotPipe& operator=(const TGnuplotPipe&) = delete;

    TGnuplotPipe(TGnuplotPipe&& other) noexcept : Pipe_(other.Pipe_) {
        other.Pipe_ = nullptr;
    }
    TGnuplotPipe& operator=(TGnuplotPipe&& other) noexcept {
        if (this != &other) {
            if (Pipe_) {
                pclose(Pipe_);
            }
            Pipe_ = other.Pipe_;
            other.Pipe_ = nullptr;
        }
        return *this;
    }

    ~TGnuplotPipe() {
        if (Pipe_) {
            pclose(Pipe_);
        }
    }

    void Write(const std::string& commands) {
        if (!Pipe_) {
            throw std::runtime_error("Pipe is closed or invalid.");
        }
        if (fwrite(commands.c_str(), sizeof(char), commands.size(), Pipe_) != commands.size()) {
            throw std::runtime_error("Failed to write to Gnuplot pipe.");
        }
    }

private:
    FILE* Pipe_;
};


inline std::string GenerateGnuplotCommands(const std::string& filename, long double left, long double right) {
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

inline void Plot(const std::string& filename, long double left, long double right) {
    std::ifstream dataFile(filename);
    if (!dataFile.is_open()) {
        throw std::runtime_error("Error: Data file '" + filename + "' not found");
    }

    try {
        const std::string commands = GenerateGnuplotCommands(filename, left, right);

        TGnuplotPipe gnuplot;
        gnuplot.Write(commands);
        #ifdef PRINT
        std::cout << "Plot successfully created: " << filename << ".png\n";
        #endif // PRINT
    } catch (const std::exception& error) {
        throw std::runtime_error(std::string("Failed to generate plot: ") + error.what());
    }
}

} // namespace NGnuplot