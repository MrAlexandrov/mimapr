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
    // TMatrix<long double> add = {
    //     {-a / L - b / 2 + c * L / 3,  a / L + b / 2 + c * L / 6}, 
    //     { a / L - b / 2 + c * L / 6, -a / L + b / 2 + c * L / 3}
    // };
    // for (int i = index; i < index + 2; ++i) {
    //     for (int j = index; j < index + 2; ++j) {
    //         resultMatrix[i][j] += add[i - index][j - index];
    //     }
    // }
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
    // TMatrix<long double> addMatrix = { 
    // {   -a *  37 / (10 * L) - b / 2       + c * 8  * L / 105,
    //          a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560,
    //         -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140,
    //          a *  13 / (40 * L) + b * 7  / 80 + c * 19 * L / 1680   },
    // {    a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560,
    //         -a *  54 / ( 5 * L)               + c * 27 * L / 70,
    //          a * 297 / (40 * L) + b * 81 / 80 - c * 27 * L / 560,
    //         -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140    },
    // {   -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140,
    //          a * 297 / (40 * L) - b * 81 / 80 - c * 27 * L / 560,
    //         -a *  54 / ( 5 * L)               + c * 27 * L / 70,
    //          a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560    },
    // {    a *  13 / (40 * L) - b * 7  / 80 + c * 19 * L / 1680,
    //         -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140,
    //          a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560,
    //         -a *  37 / (10 * L) + b / 2       + c * 8  * L / 105    }
    // };
    // for (int i = index; i < index + 4; ++i) {
    //     for (int j = index; j < index + 4; ++j) {
    //         resultMatrix[i][j] += addMatrix[i - index][j - index];
    //     }
    // }

    
    // TMatrix<> aCoef = {
    //     {   -37.0   /   10.0,   189.0   /   40.0,   27.0    /   20.0,   13.0    /   40.0},
    //     {   189.0   /   40.0,   -54.0   /   5.0,    297.0   /   40.0,   -27.0   /   20.0},
    //     {   -27.0   /   20.0,   297.0   /   40.0,   -54.0   /   5.0,    189.0   /   40.0},
    //     {   13.0    /   40.0,   -27.0   /   20.0,   189.0   /   40.0,   -37.0   /   10.0}
    // };

    // TMatrix<> bCoef = {
    //     {   -1.0    /   2.0,    57.0    /   80.0,   -3.0    /   10.0,   7.0     /   80.0},
    //     {   -57.0   /   80.0,   0.0             ,   81.0    /   80.0,   -3.0    /   10.0},
    //     {   3.0     /   10.0,   -81.0   /   80.0,   0.0             ,   57.0    /   80.0},
    //     {   -7.0    /   80.0,   3.0     /   10.0,   -57.0   /   80.0,   1.0     /   2.0}
    // };

    // TMatrix<> cCoef = {
    //     {   8.0     /   105.0,  33.0    /   560.0,  -3.0    /   140.0,  19.0    /   1680.0},
    //     {   33.0    /   27.0,   -27.0   /   560.0,  -3.0    /   140.0,  -3.0    /   140.0},
    //     {   -3.0    /   140.0,  -27.0   /   560.0,  27.0    /   70.0,   33.0    /   560.0},
    //     {   19.0    /   1680.0, -3.0    /   140.0,  33.0    /   560.0,  8.0     /   105.0}
    // };

    // for (int i = 0; i < 4; ++i) {
    //     for (int j = 0; j < 4; ++j) {
    //         resultMatrix[index + i][index + j] += 
    //             a * aCoef[i][j] * L + 
    //             b * bCoef[i][j]     +
    //             c * cCoef[i][j] * L;
    //     }
    // }

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

void applyRestrictions(TMatrix<long double>& matrix, 
                       TMatrix<long double>& vector, 
                       const Restriction& restriction, 
                       const int row, const long double coef) {
    switch (restriction.grade) {
        case RestrictGrade::First: {
            for (int col = 0; col < matrix.cols(); ++col) {
                matrix[row][col] = 0;
            }
            matrix[row][row] = 1;
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

void solveSLAU(TMatrix<long double>& xVector, TMatrix<long double>& aMatrix, TMatrix<long double>& bVector) {
    int size = aMatrix.rows();
    for (int i = 0; i < size - 1; ++i) {
        for (int j = i + 1; j < size; ++j) {
            long double factor = aMatrix[j][i] / aMatrix[i][i];
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

} // namespace

int main(int argc, char* argv[]) {
    auto opt = getOptions(argc, argv);
    if (!opt.has_value() || opt->help) {
        return 0;
    }

    constexpr long double a = 1, b = 0, c = -1, d = -10;
    constexpr Restriction lower = {RestrictGrade::Second, 2, 10};
    constexpr Restriction upper = {RestrictGrade::First, 8, 5};

    int size = opt->elemAmount * (opt->type == ElementType::Linear ? 1 : 3) + 1;
    long double step = static_cast<long double>(upper.pos - lower.pos) / opt->elemAmount;

    TMatrix<long double> stiffnessMatrix(size, size, 0.0);
    TMatrix<long double> loadVector(size, 1, 0.0);
    TMatrix<long double> displacements(size, 1, 0.0);

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
    
    long double maxError = 0.0;
    std::vector<long double> nodes(size);
    std::vector<long double> errors(size);
    std::vector<long double> displacementsReal(size);

    long double nodePosition = lower.pos;
    for (int i = 0; i < size; ++i) {
        long double realValue = realSolve(nodePosition);
        long double error = std::fabs(displacements[i][0] - realValue);
        maxError = std::fmax(maxError, error);

        nodes[i] = nodePosition;
        errors[i] = error;
        displacementsReal[i] = realValue;

        nodePosition += step / (opt->type == ElementType::Linear ? 1 : 3);
    }

    // #ifdef PRINT
    std::cout << "\nMax error: " << maxError << "\n";
    // #endif // PRINT
    try {
        std::string filename = (opt->type == ElementType::Linear ? "linear" : "cubic") 
                       + std::string("_") 
                       + std::to_string(opt->elemAmount);

        std::ofstream outFile(filename + ".txt");
        if (!outFile.is_open()) {
            throw std::ios_base::failure("Failed to open the output file: " + filename + ".txt");
        }
        constexpr long double EPS = 1e-12;
        for (int i = 0; i < size; ++i) {
            outFile << std::setw(4)  << nodes[i]
                    << std::setw(12) << (-EPS < displacementsReal[i] && displacementsReal[i] < EPS ? 0 : displacementsReal[i])
                    << std::setw(12) << (-EPS < displacements[i][0] && displacements[i][0] < EPS ? 0 : displacements[i][0])
                    << std::setw(12) << (-EPS < errors[i] && errors[i] < EPS ? 0 : errors[i]) << '\n';
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
