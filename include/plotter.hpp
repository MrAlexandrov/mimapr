#pragma once
#include <cassert>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

namespace NPlotter {

namespace {

class TGraph {
public:
    TGraph(const std::string& title, const std::vector<long double> data)
        : Title_(title)
        , Data_(data)
        {}

    void SetTitle(std::string&& title) {
        Title_ = std::move(title);
    }

    void SetData(const std::vector<long double>& data) {
        Data_ = data;
    }

    std::string GetTitle() const {
        return Title_;
    }

    std::vector<long double> GetData() const {
        return Data_;
    }

    long double GetData(int index) const {
        assert(0 <= index && index < Data_.size());
        return Data_[index];
    }

private:
    std::string Title_;
    std::vector<long double> Data_;
};

} // namespace

class TPlotter final {
public:
    explicit TPlotter(const std::string& filename)
        : ImageName_(filename + ".png")
        , DataFile_(filename + ".csv")
        , Width_(800)
        , Height_(600)
        , Title_("Displacements and Errors")
        , XLabel_("Node position")
        , YLabel_("Values")
        , XRangeLeft_(0.0)
        , XRangeRight_(10.0)
        {}

    void SetXValues(const std::vector<long double> values) {
        assert(Datas_.empty() || Datas_.front().GetData().size() == values.size());
        XValues_ = values;
    }

    void SetXRangeLeft(long double value) {
        XRangeLeft_ = value;
    }

    void SetXRangeRight(long double value) {
        XRangeRight_ = value;
    }

    void AddGraphic(std::string&& name, std::vector<long double> data) {
        assert(XValues_.empty() || XValues_.size() == data.size());
        assert(Datas_.empty() || Datas_.front().GetData().size() == data.size());
        Datas_.emplace_back(name, data);
    }

    void Plot() const {
        this->SaveDataToFile();
        try {
            FILE* pipe{popen("gnuplot -persist", "w")};
            if (!pipe) {
                throw std::runtime_error("Failed to open pipe to gnuplot.");
            }
            std::string commands = GenerateGnuplotCommands();
            if (fwrite(commands.c_str(), sizeof(char), commands.size(), pipe) != commands.size()) {
                throw std::runtime_error("Failed to write to Gnuplot pipe.");
            }
            pclose(pipe);
        } catch (std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknow error, while plotting" << std::endl;
        }
    }

private:
    long double PrintValue(long double value) const {
        constexpr long double EPS = 1e-12; 
        return (-EPS < value && value < EPS ? 0 : value);
    }

    void SaveDataToFile() const {
        int size = Datas_.front().GetData().size();
        std::ofstream outFile(DataFile_);
        if (!outFile.is_open()) {
            throw std::ios_base::failure("Failed to open the output file: " + DataFile_);
        }
        for (int row = 0; row < size; ++row) {
            outFile << std::setw(4) << PrintValue(XValues_[row]);
            for (auto& data : Datas_) {
                outFile << std::setw(12) << PrintValue(data.GetData(row));
            }
            outFile << "\n";
        }
        outFile.flush();
        outFile.close();
    }

private:
    std::string GenerateGnuplotCommands() const {
        std::ostringstream commands;
        commands << "set term png size " << Width_ << "," << Height_ << "\n";
        commands << "set output '" << ImageName_ << "'\n";
        commands << "set title '" << Title_ << "'\n";
        commands << "set xlabel '" << XLabel_ << "'\n";
        commands << "set ylabel '" << YLabel_ << "'\n";
        commands << "set xrange [" << XRangeLeft_ << ":" << XRangeRight_ << "]\n";
        commands << "plot ";

        int index = 2;
        for (const auto& data : Datas_) {
            commands << "'" << DataFile_ << "' using 1:" << index++ << " with lines title '" << data.GetTitle() << "', ";
        }
        std::string result = commands.str();
        result.pop_back();
        result.pop_back();
        return result;
    }

private:
    unsigned int Width_;
    unsigned int Height_;
    std::string ImageName_;
    std::string DataFile_;
    std::string Title_;
    std::string XLabel_;
    std::string YLabel_;
    long double XRangeLeft_;
    long double XRangeRight_;
    std::vector <long double> XValues_;
    std::vector<TGraph> Datas_;
};

} // NPlotter