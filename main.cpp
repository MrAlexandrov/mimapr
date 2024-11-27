#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <cassert>

#include "initialize.hpp"
#include "matrix.hpp"
#include "gnuplot.hpp"

namespace {

using namespace NMatrix;
using namespace NGnuplot;
using namespace NOptions;

void MakeLinearSLAU(TMatrix<>& resultMatrix, 
                    TMatrix<>& resultVector, 
                    const long double L,
                    const int index)
{
    resultMatrix[index][index]          += -A / L   - B / 2 + C * L / 3;
    resultMatrix[index][index + 1]      +=  A / L   + B / 2 + C * L / 6;
    resultMatrix[index + 1][index]      +=  A / L   - B / 2 + C * L / 6;
    resultMatrix[index + 1][index + 1]  += -A / L   + B / 2 + C * L / 3;

    resultVector[index][0]              += -D * L   / 2;
    resultVector[index + 1][0]          += -D * L   / 2;
}

void MakeCubicSLAU(TMatrix<>& resultMatrix, 
                   TMatrix<>& resultVector, 
                   const long double L,
                   const int index) 
{   
    resultMatrix[index][index]          += -A *  37 / (10 * L) - B / 2       + C * 8  * L / 105;
    resultMatrix[index][index + 1]      +=  A * 189 / (40 * L) + B * 57 / 80 + C * 33 * L / 560;
    resultMatrix[index][index + 2]      += -A *  27 / (20 * L) - B * 3  / 10 - C * 3  * L / 140;
    resultMatrix[index][index + 3]      +=  A *  13 / (40 * L) + B * 7  / 80 + C * 19 * L / 1680;

    resultMatrix[index + 1][index]      +=  A * 189 / (40 * L) - B * 57 / 80 + C * 33 * L / 560;
    resultMatrix[index + 1][index + 1]  += -A *  54 / (5  * L)               + C * 27 * L / 70;
    resultMatrix[index + 1][index + 2]  +=  A * 297 / (40 * L) + B * 81 / 80 - C * 27 * L / 560;
    resultMatrix[index + 1][index + 3]  += -A *  27 / (20 * L) - B * 3  / 10 - C * 3  * L / 140;

    resultMatrix[index + 2][index]      += -A *  27 / (20 * L) + B * 3  / 10 - C * 3  * L / 140;
    resultMatrix[index + 2][index + 1]  +=  A * 297 / (40 * L) - B * 81 / 80 - C * 27 * L / 560;
    resultMatrix[index + 2][index + 2]  += -A *  54 / (5  * L)               + C * 27 * L / 70;
    resultMatrix[index + 2][index + 3]  +=  A * 189 / (40 * L) + B * 57 / 80 + C * 33 * L / 560;

    resultMatrix[index + 3][index]      +=  A *  13 / (40 * L) - B * 7  / 80 + C * 19 * L / 1680;
    resultMatrix[index + 3][index + 1]  += -A *  27 / (20 * L) + B * 3  / 10 - C * 3  * L / 140;
    resultMatrix[index + 3][index + 2]  +=  A * 189 / (40 * L) - B * 57 / 80 + C * 33 * L / 560;
    resultMatrix[index + 3][index + 3]  += -A *  37 / (10 * L) + B / 2       + C * 8  * L / 105;

    resultVector[index][0]              += -D * L / 8;
    resultVector[index + 1][0]          += -D * 3 * L / 8;
    resultVector[index + 2][0]          += -D * 3 * L / 8;
    resultVector[index + 3][0]          += -D * L / 8;
}

void ApplyRestriction(TMatrix<>& matrix, 
                      TMatrix<>& vector, 
                      const TRestriction& restriction, 
                      const int row, const long double coef) {
    switch (restriction.Grade) {
        case ERestrictionGrade::FIRST: {
            for (int col = 0, end = matrix.cols(); col < end; ++col) {
                matrix[row][col] = 0;
            }
            matrix[row][row] = 1;
            vector[row][0] = restriction.Value;
            break;
        }
        case ERestrictionGrade::SECOND: {
            vector[row][0] += coef * (row == matrix.rows() - 1 ? -restriction.Value : restriction.Value);
            break;
        }
        case ERestrictionGrade::THIRD: {
            matrix[row][row] += coef * restriction.Value;
            break;
        }
        default: {
            throw std::invalid_argument("Unknown restriction grade.");
        }
    }
}

void SolveSLAU(TMatrix<>& xVector, TMatrix<>& aMatrix, TMatrix<>& bVector) {
    int size = aMatrix.rows();
    for (int i = 0; i < size - 1; ++i) {
        long double baseElement = aMatrix[i][i];
        for (int j = i + 1; j < size; ++j) {
            long double factor = aMatrix[j][i] / baseElement;
            for (int k = i; k < size; ++k) {
                aMatrix[j][k] -= factor * aMatrix[i][k];
            }
            bVector[j][0] -= factor * bVector[i][0];
        }
    }

    for (int i = size - 1; i >= 0; --i) {
        long double sum = 0;
        for (int j = i + 1; j < size; ++j) {
            sum += aMatrix[i][j] * xVector[j][0];
        }
        xVector[i][0] = (bVector[i][0] - sum) / aMatrix[i][i];
    }
}

long double RealSolve(const long double x) {
	long double C1 = -(5 * expl(8) * (-3 + 2 * expl(6))) / (1 + expl(12));
	long double C2 = 5 * (2 + 3 * expl(6)) / (expl(2) * (1 + expl(12)));
	return expl(-x) * C1 + expl(x) * C2 - 10;
}

void SaveResultsToFile(const std::string& filename, 
                       const std::vector<long double>& nodes,
                       const std::vector<long double>& displacementsReal,
                       const TMatrix<>&     displacements,
                       const std::vector<long double>& errors) 
{
    try {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            throw std::ios_base::failure("Failed to open the output file: " + filename);
        }
        constexpr long double EPS = 1e-12;
        for (size_t i = 0, end = nodes.size(); i < end; ++i) {
            outFile << std::setw(4) << nodes[i]
                    << std::setw(12) << (-EPS < displacementsReal[i] && displacementsReal[i] < EPS ? 0 : displacementsReal[i])
                    << std::setw(12) << (-EPS < displacements[i][0] && displacements[i][0] < EPS ? 0 : displacements[i][0])
                    << std::setw(12) << (-EPS < errors[i] && errors[i] < EPS ? 0 : errors[i]) << '\n';
        }
        outFile.close();
        std::cout << "Results successfully saved to " << filename << "\n";
    } catch (const std::exception& error) {
        std::cerr << "An error occurred while saving results: " << error.what() << "\n";
    }
}

} // namespace

int main(int argc, char* argv[]) {
    const std::optional<TOptions> opt = GetOptions(argc, argv);
    if (!opt.has_value() || opt->Help) {
        return 0;
    }

    const std::vector<TRestriction> restrictions = {
        LOWER,
        UPPER,
    };

    const long double minimum = std::min_element(restrictions.begin(), restrictions.end())->Position;
    const long double maximum = std::max_element(restrictions.begin(), restrictions.end())->Position;

    const int size = opt->ElementsAmount * (opt->Type == EElementType::LINEAR ? 1 : 3) + 1;
    const long double step = static_cast<long double>(maximum - minimum) / opt->ElementsAmount;

    TMatrix<> stiffnessMatrix(size, size, 0.0);
    TMatrix<> loadVector(size, 1, 0.0);
    TMatrix<> displacements(size, 1, 0.0);

    switch (opt->Type) {
        case NOptions::EElementType::LINEAR: {
            for (int i = 0; i < size - 1; ++i) {
                MakeLinearSLAU(stiffnessMatrix, loadVector, step, i);
            }
            break;
        }
        case NOptions::EElementType::CUBIC: {
            for (int i = 0; i < size - 3; i += 3) {
                MakeCubicSLAU(stiffnessMatrix, loadVector, step, i);
            }
            break;
        }
        default: {
            throw std::runtime_error("Unknown element type");
        }
    }

    auto getRowByPosition = [&](long double position) -> int {
        assert(minimum <= position && position <= maximum);
        long double shift = position - minimum;
        if (opt->Type == EElementType::CUBIC) {
            shift *= 3;
        }
        int row = std::round(shift / step);
        row = std::max(0, row);
        row = std::min(size - 1, row);
        return row;
    };

    for (auto&& current : restrictions) {
        ApplyRestriction(stiffnessMatrix, loadVector, current, getRowByPosition(current.Position), A);
    }

    SolveSLAU(displacements, stiffnessMatrix, loadVector);
    
    #ifdef PRINT
    for (int i = 0; i < size; ++i) {
        std::cout << "U(" << i << ") = " << displacements[i][0] << "\n";
    }
    #endif // PRINT
    
    long double maxError = 0.0;
    std::vector<long double> nodes(size);
    std::vector<long double> errors(size);
    std::vector<long double> displacementsReal(size);

    long double nodePosition = minimum;
    for (int i = 0; i < size; ++i) {
        long double realValue = RealSolve(nodePosition);
        long double error = std::fabs(displacements[i][0] - realValue);
        maxError = std::fmax(maxError, error);

        nodes[i] = nodePosition;
        errors[i] = error;
        displacementsReal[i] = realValue;

        nodePosition += step / (opt->Type == EElementType::LINEAR ? 1 : 3);
    }

    // #ifdef PRINT
    std::cout << "\nMax error: " << maxError << "\n";
    // #endif // PRINT

    const std::string filename = (opt->Type == EElementType::LINEAR ? "linear" : "cubic") 
                    + std::string("_") 
                    + std::to_string(opt->ElementsAmount)
                    + std::string(".txt");

    SaveResultsToFile(filename, nodes, displacementsReal, displacements, errors);
    Plot(filename, minimum, maximum);

    return 0;
}
