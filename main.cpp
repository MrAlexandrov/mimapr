#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <cassert>

#include "matrix.hpp"
#include "options.hpp"
#include "gnuplot.hpp"

namespace {

using namespace NMatrix;
using namespace NGnuplot;
using namespace NOptions;

void makeLinearSLAU(TMatrix<long double>& resultMatrix, 
                    TMatrix<long double>& resultVector, 
                    const long double L, 
                    const long double a, const long double b, const long double c, const long double d, 
                    const int index) 
{
    resultMatrix[index][index]          += -a / L   - b / 2 + c * L / 3;
    resultMatrix[index][index + 1]      +=  a / L   + b / 2 + c * L / 6;
    resultMatrix[index + 1][index]      +=  a / L   - b / 2 + c * L / 6;
    resultMatrix[index + 1][index + 1]  += -a / L   + b / 2 + c * L / 3;

    resultVector[index][0]          += -d * L   / 2;
    resultVector[index + 1][0]      += -d * L   / 2;
}

void makeCubicSLAU(TMatrix<long double>& resultMatrix, 
                   TMatrix<long double>& resultVector, 
                   const long double L, 
                   const long double a, const long double b, const long double c, const long double d, 
                   const int index) 
{   
    resultMatrix[index][index]          += -a *  37 / (10 * L) - b / 2       + c * 8  * L / 105;
    resultMatrix[index][index + 1]      +=  a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[index][index + 2]      += -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[index][index + 3]      +=  a *  13 / (40 * L) + b * 7  / 80 + c * 19 * L / 1680;

    resultMatrix[index + 1][index]      +=  a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[index + 1][index + 1]  += -a *  54 / (5  * L)               + c * 27 * L / 70;
    resultMatrix[index + 1][index + 2]  +=  a * 297 / (40 * L) + b * 81 / 80 - c * 27 * L / 560;
    resultMatrix[index + 1][index + 3]  += -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140;

    resultMatrix[index + 2][index]      += -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[index + 2][index + 1]  +=  a * 297 / (40 * L) - b * 81 / 80 - c * 27 * L / 560;
    resultMatrix[index + 2][index + 2]  += -a *  54 / (5  * L)               + c * 27 * L / 70;
    resultMatrix[index + 2][index + 3]  +=  a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560;

    resultMatrix[index + 3][index]      +=  a *  13 / (40 * L) - b * 7  / 80 + c * 19 * L / 1680;
    resultMatrix[index + 3][index + 1]  += -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[index + 3][index + 2]  +=  a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[index + 3][index + 3]  += -a *  37 / (10 * L) + b / 2       + c * 8  * L / 105;

    resultVector[index][0]          += -d * L / 8;
    resultVector[index + 1][0]      += -d * 3 * L / 8;
    resultVector[index + 2][0]      += -d * 3 * L / 8;
    resultVector[index + 3][0]      += -d * L / 8;
}

void applyRestriction(TMatrix<long double>& matrix, 
                       TMatrix<long double>& vector, 
                       const TRestriction& restriction, 
                       const int row, const long double coef) {
    switch (restriction.grade) {
        case ERestrictionGrade::First: {
            for (int col = 0, end = matrix.cols(); col < end; ++col) {
                matrix[row][col] = 0;
            }
            matrix[row][row] = 1;
            vector[row][0] = restriction.value;
            break;
        }
        case ERestrictionGrade::Second: {
            vector[row][0] += coef * (row == 0 ? restriction.value : (row == matrix.rows() - 1 ? -restriction.value : restriction.value));
            break;
        }
        case ERestrictionGrade::Third: {
            matrix[row][row] += coef * restriction.value;
            break;
        }
        default: {
            throw std::invalid_argument("Unknown restriction grade.");
        }
    }
}

void solveSLAU(TMatrix<long double>& xVector, TMatrix<long double>& aMatrix, TMatrix<long double>& bVector) {
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

long double realSolve(const long double x) {
	long double C1 = -(5 * expl(8) * (-3 + 2 * expl(6))) / (1 + expl(12));
	long double C2 = 5 * (2 + 3 * expl(6)) / (expl(2) * (1 + expl(12)));
	return expl(-x) * C1 + expl(x) * C2 - 10;
}

void SaveResultsToFile(const std::string& filename, 
                       const std::vector<long double>& nodes,
                       const std::vector<long double>& displacementsReal,
                       const TMatrix<long double>&     displacements,
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
    auto opt = getOptions(argc, argv);
    if (!opt.has_value() || opt->help) {
        return 0;
    }

    constexpr long double a = 1, b = 0, c = -1, d = -10;
    constexpr TRestriction lower = {ERestrictionGrade::Second, 2, 10};
    constexpr TRestriction upper = {ERestrictionGrade::First, 8, 5};

    std::vector<TRestriction> Restrictions = {
        lower,
        upper,
    };

    int minimum = (*std::min_element(Restrictions.begin(), Restrictions.end())).position;
    int maximum = (*std::max_element(Restrictions.begin(), Restrictions.end())).position;

    int size = opt->elementsAmount * (opt->type == EElementType::Linear ? 1 : 3) + 1;
    long double step = static_cast<long double>(maximum - minimum) / opt->elementsAmount;

    TMatrix<long double> stiffnessMatrix(size, size, 0.0);
    TMatrix<long double> loadVector(size, 1, 0.0);
    TMatrix<long double> displacements(size, 1, 0.0);

    switch (opt->type) {
        case NOptions::EElementType::Linear: {
            for (int i = 0; i < size - 1; ++i) {
                makeLinearSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
            }
            break;
        }
        case NOptions::EElementType::Cubic: {
            for (int i = 0; i < size - 3; i += 3) {
                makeCubicSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
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
        if (opt->type == EElementType::Cubic) {
            shift *= 3;
        }
        int row = std::round(shift / step);
        row = std::max(0, row);
        row = std::min(size - 1, row);
        return row;
    };

    for (auto&& current : Restrictions) {
        applyRestriction(stiffnessMatrix, loadVector, current, getRowByPosition(current.position), a);
    }

    solveSLAU(displacements, stiffnessMatrix, loadVector);
    
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
        long double realValue = realSolve(nodePosition);
        long double error = std::fabs(displacements[i][0] - realValue);
        maxError = std::fmax(maxError, error);

        nodes[i] = nodePosition;
        errors[i] = error;
        displacementsReal[i] = realValue;

        nodePosition += step / (opt->type == EElementType::Linear ? 1 : 3);
    }

    // #ifdef PRINT
    std::cout << "\nMax error: " << maxError << "\n";
    // #endif // PRINT

    std::string filename = (opt->type == EElementType::Linear ? "linear" : "cubic") 
                    + std::string("_") 
                    + std::to_string(opt->elementsAmount)
                    + std::string(".txt");

    SaveResultsToFile(filename, nodes, displacementsReal, displacements, errors);
    Plot(filename, minimum, maximum);

    return 0;
}
