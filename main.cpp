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

void makeLinearSLAU(TMatrix<>& resultMatrix, 
                    TMatrix<>& resultVector, 
                    double L, 
                    double a, double b, double c, double d, 
                    int i) 
{
    resultMatrix[i][i]          += -a / L   - b / 2 + c * L / 3;
    resultMatrix[i][i + 1]      +=  a / L   + b / 2 + c * L / 6;
    resultMatrix[i + 1][i]      +=  a / L   - b / 2 + c * L / 6;
    resultMatrix[i + 1][i + 1]  += -a / L   + b / 2 + c * L / 3;

    resultVector[i][0]          += -d * L   / 2;
    resultVector[i + 1][0]      += -d * L   / 2;
}

void makeCubicSLAU(TMatrix<>& resultMatrix, 
                   TMatrix<>& resultVector, 
                   double L, 
                   double a, double b, double c, double d, 
                   int i) 
{
    resultMatrix[i][i]          += -a *  37 / (10 * L) - b / 2       + c * 8  * L / 105;
    resultMatrix[i][i + 1]      +=  a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[i][i + 2]      += -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[i][i + 3]      +=  a *  13 / (40 * L) + b * 7  / 80 + c * 19 * L / 1680;

    resultMatrix[i + 1][i]      +=  a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[i + 1][i + 1]  += -a *  54 / (5  * L)               + c * 27 * L / 70;
    resultMatrix[i + 1][i + 2]  +=  a * 297 / (40 * L) + b * 81 / 80 - c * 27 * L / 560;
    resultMatrix[i + 1][i + 3]  += -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140;

    resultMatrix[i + 2][i]      += -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[i + 2][i + 1]  +=  a * 297 / (40 * L) - b * 81 / 80 - c * 27 * L / 560;
    resultMatrix[i + 2][i + 2]  += -a *  54 / (5  * L)               + c * 27 * L / 70;
    resultMatrix[i + 2][i + 3]  +=  a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560;

    resultMatrix[i + 3][i]      +=  a *  13 / (40 * L) - b * 7  / 80 + c * 19 * L / 1680;
    resultMatrix[i + 3][i + 1]  += -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[i + 3][i + 2]  +=  a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[i + 3][i + 3]  += -a *  37 / (10 * L) + b / 2       + c * 8  * L / 105;

    resultVector[i][0]          += -d * L / 8;
    resultVector[i + 1][0]      += -d * 3 * L / 8;
    resultVector[i + 2][0]      += -d * 3 * L / 8;
    resultVector[i + 3][0]      += -d * L / 8;
}

void applyRestrictions(TMatrix<>& matrix, TMatrix<>& vector, Restriction restriction, int row, double coef) {
    switch (restriction.grade) {
        case RestrictGrade::First: {
            for (int col = 0; col < matrix.cols(); ++col) {
                matrix[row][col] = 0.0;
            }
            matrix[row][row] = 1.0;
            vector[row][0] = restriction.val;
            break;
        }
        case RestrictGrade::Second: {
            vector[row][0] += coef * (row == 0 ? restriction.val : (row == matrix.rows() - 1 ? -restriction.val : 0));
            break;
        }
        case RestrictGrade::Third: {
            matrix[row][row] += coef * restriction.val;
            break;
        }
        default:
            throw std::invalid_argument("Unknown restriction grade.");
    }
}

void solveSLAU(TMatrix<>& xVector, TMatrix<>& aMatrix, TMatrix<>& bVector) {
    int size = aMatrix.rows();
    for (int i = 0; i < size - 1; ++i) {
        for (int j = i + 1; j < size; ++j) {
            double factor = aMatrix[j][i] / aMatrix[i][i];
            for (int k = i; k < size; ++k) {
                aMatrix[j][k] -= factor * aMatrix[i][k];
            }
            bVector[j][0] -= factor * bVector[i][0];
        }
    }

    for (int i = size - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < size; ++j) {
            sum += aMatrix[i][j] * xVector[j][0];
        }
        xVector[i][0] = (bVector[i][0] - sum) / aMatrix[i][i];
    }
}

double realSolve(double x) {
	double C1 = -(5 * exp(8) * (-3 + 2 * exp(6))) / (1 + exp(12));
	double C2 = 5 * (2 + 3 * exp(6)) / (exp(2) * (1 + exp(12)));
	return exp(-x) * C1 + exp(x) * C2 - 10;
}

} // namespace

int main(int argc, char* argv[]) {
    auto opt = getOptions(argc, argv);
    if (!opt.has_value() || opt->help) {
        return 0;
    }

    const double a = 1, b = 0, c = -1, d = -10;
    const Restriction lower = {RestrictGrade::Second, 2, 10};
    const Restriction upper = {RestrictGrade::First, 8, 5};

    int size = opt->elemAmount * (opt->type == ElementType::Linear ? 1 : 3) + 1;
    double step = (upper.pos - lower.pos) / opt->elemAmount;

    TMatrix<> stiffnessMatrix(size, size, 0.0);
    TMatrix<> loadVector(size, 1, 0.0);
    TMatrix<> displacements(size, 1, 0.0);

    if (opt->type == ElementType::Linear) {
        for (int i = 0; i < size - 1; ++i) {
            makeLinearSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
        }
    } else if (opt->type == ElementType::Cubic) {
        for (int i = 0; i < size - 3; i += 3) {
            makeCubicSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
        }
    } else {
        throw std::runtime_error("ElementType::Unknown");
    }
    applyRestrictions(stiffnessMatrix, loadVector, lower, 0, a);
    applyRestrictions(stiffnessMatrix, loadVector, upper, size - 1, a);
    solveSLAU(displacements, stiffnessMatrix, loadVector);
    
    #ifdef PRINT
    for (int i = 0; i < size; ++i) {
        std::cout << "U(" << i << ") = " << displacements[i][0] << "\n";
    }
    #endif // PRINT
    
    double maxError = 0.0;
    std::vector<double> nodes(size);
    std::vector<double> errors(size);
    std::vector<double> displacementsReal(size);

    double nodePosition = lower.pos;
    for (int i = 0; i < size; ++i) {
        double realValue = realSolve(nodePosition);
        double error = std::fabs(displacements[i][0] - realValue);
        maxError = std::max(maxError, error);

        nodes[i] = nodePosition;
        errors[i] = error;
        displacementsReal[i] = realValue;

        nodePosition += step / (opt->type == ElementType::Linear ? 1 : 3);
    }

    #ifdef PRINT
    std::cout << "\nMax error: " << maxError << "\n";
    #endif // PRINT
    try {
        std::string filename = (opt->type == ElementType::Linear ? "linear" : "cubic") 
                       + std::string("_") 
                       + std::to_string(opt->elemAmount);

        std::ofstream outFile(filename + ".txt");
        if (!outFile) {
            throw std::ios_base::failure("Failed to open the output file: " + filename + ".txt");
        }

        for (int i = 0; i < size; ++i) {
            outFile << std::setw(10) << nodes[i]
                    << std::setw(10) << displacementsReal[i]
                    << std::setw(10) << displacements[i][0]
                    << std::setw(10) << errors[i] << '\n';
        }
        outFile.flush();
        #ifdef PRINT
        std::cout << "Results successfully saved to " << filename << ".txt\n";
        #endif // PRINT
        plot(filename, lower.pos, upper.pos);
    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << '\n';
    }

    return 0;
}
