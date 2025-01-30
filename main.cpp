#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <filesystem>

#include "initialize.hpp"
#include "matrix.hpp"
#include "options.hpp"
#include "plotter.hpp"

namespace {

using namespace NMatrix;
using namespace NPlotter;
using namespace NOptions;

void MakeLinearSLAU(
    TMatrix<double>&            resultMatrix,
    TMatrix<double>&            resultVector,
    const double                L,
    const int                   index
) {
    resultMatrix[index][index]          += -A / L   - B / 2 + C * L / 3;
    resultMatrix[index][index + 1]      +=  A / L   + B / 2 + C * L / 6;
    resultMatrix[index + 1][index]      +=  A / L   - B / 2 + C * L / 6;
    resultMatrix[index + 1][index + 1]  += -A / L   + B / 2 + C * L / 3;

    resultVector[index][0]              += -D * L   / 2;
    resultVector[index + 1][0]          += -D * L   / 2;
}

void MakeCubicSLAU(
    TMatrix<double>&            resultMatrix,
    TMatrix<double>&            resultVector,
    const double                L,
    const int                   index
) {
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

void InitializeSLAU(
    TMatrix<double>&            stiffnessMatrix,
    TMatrix<double>&            loadVector,
    const EElementType&         type,
    const int                   size,
    const double                step
) {
    switch (type) {
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
}

void ApplyRestriction(
    TMatrix<double>&          matrix,
    TMatrix<double>&          vector,
    const TRestriction&       restriction,
    const int                 row
) {
    assert(row == 0 || row == matrix.Rows() - 1);
    switch (restriction.Grade) {
        case ERestrictionGrade::FIRST: {
            for (int col = 0, end = matrix.Cols(); col < end; ++col) {
                matrix[row][col] = 0;
            }
            matrix[row][row] = 1;
            vector[row][0] = restriction.Value;
            break;
        }
        case ERestrictionGrade::SECOND: {
            vector[row][0] += A * restriction.Value;
            break;
        }
        case ERestrictionGrade::THIRD: {
            matrix[row][row] += A * restriction.Value;
            break;
        }
        default: {
            throw std::invalid_argument("Unknown restriction grade.");
        }
    }
}

void SolveSLAU(TMatrix<double>& xVector, TMatrix<double>& aMatrix, TMatrix<double>& bVector) {
    TMatrix<double> inverseMatrix = NMatrix::InvertMatrix(aMatrix);
    
    xVector = inverseMatrix * bVector;
}

// void SolveSLAU(TMatrix<double>& xVector, TMatrix<double>& aMatrix, TMatrix<double>& bVector) {
//     int size = aMatrix.Rows();
//     for (int i = 0; i < size - 1; ++i) {
//         double baseElement = aMatrix[i][i];
//         for (int j = i + 1; j < size; ++j) {
//             double factor = aMatrix[j][i] / baseElement;
//             for (int k = i; k < size; ++k) {
//                 aMatrix[j][k] -= factor * aMatrix[i][k];
//             }
//             bVector[j][0] -= factor * bVector[i][0];
//         }
//     }

//     for (int i = size - 1; i >= 0; --i) {
//         double sum = 0;
//         for (int j = i + 1; j < size; ++j) {
//             sum += aMatrix[i][j] * xVector[j][0];
//         }
//         xVector[i][0] = (bVector[i][0] - sum) / aMatrix[i][i];
//     }
// }

inline double RealSolve(const double x) {
    double C1 = -(5 * expl(8) * (-3 + 2 * expl(6))) / (1 + expl(12));
    double C2 = 5 * (2 + 3 * expl(6)) / (expl(2) * (1 + expl(12)));
    return expl(-x) * C1 + expl(x) * C2 - 10;
}

inline double RealSolveAdditional(const double x) {
    return static_cast<double>(5) * (expl(8 - x) + 2 * expl(x - 2) + expl(x + 4) - 2);
}

std::vector<double> FillNodes(
    double minimum,
    double size,
    double add
) {
    std::vector<double> result;
    double current = minimum;
    for (int i = 0; i < size; ++i, current += add) {
        result.push_back(current);
    }
    return result;
}

std::vector<double> FillRealValues(
    const std::vector<double>& positions, 
    const auto& function
) {
    std::vector<double> result;
    for (const auto& x : positions) {
        result.emplace_back(function(x));
    }
    return result;
}

std::vector<double> FillErrors(
    const std::vector<double>& solution,
    const std::vector<double>& realValues
) {
    assert(solution.size() == realValues.size());
    std::vector<double> result;
    for (int i = 0, end = solution.size(); i < end; ++i) {
        result.emplace_back(std::fabs(solution[i] - realValues[i]));
    }
    return result;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& other) {
    for (const auto& i : other) {
        out << i << ' ';
    }
    return out;
}

namespace fs = std::filesystem;

std::vector<int> sizes = {/*1, 3, 20, 40, 100, 200, 400, */800, 1600, 3200};
std::vector<std::string> ke_types = {"linear", "cubic"};

void MoveFilesToFolder(const std::string& extension, const std::string& folder) {
    fs::create_directories(folder);
    for (const auto& entry : fs::directory_iterator(".")) {
        if (entry.path().extension() == extension) {
            fs::rename(entry.path(), folder + "/" + entry.path().filename().string());
        }
    }
}

} // anonimus namespace

int main() {
    for (const auto& size : sizes) {
        for (const auto& ke_type : ke_types) {
            if ((ke_type == "linear" && size == 1) || (ke_type == "cubic" && size == 3)) {
                continue;
            }

            TOptions options;
            options.Type = (ke_type == "linear") ? EElementType::LINEAR : EElementType::CUBIC;
            options.ElementsAmount = size;

            const int matrixSize = options.ElementsAmount * EElementTypeToInt[options.Type] + 1;
            const double step = static_cast<double>(maximum - minimum) / options.ElementsAmount;

            TMatrix<double> stiffnessMatrix(matrixSize, matrixSize, 0.0);
            TMatrix<double> loadVector(matrixSize, 1, 0.0);
            TMatrix<double> displacements(matrixSize, 1, 0.0);

            InitializeSLAU(stiffnessMatrix, loadVector, options.Type, matrixSize, step);

            auto getRowByPosition = [&](double position) -> int {
                assert(minimum <= position && position <= maximum);
                double shift = position - minimum;
                shift *= EElementTypeToInt[options.Type];
                int row = std::round(shift / step);
                return std::clamp(row, 0, matrixSize - 1);
            };

            for (const auto& current : restrictions) {
                ApplyRestriction(stiffnessMatrix, loadVector, current, getRowByPosition(current.Position));
            }

            SolveSLAU(displacements, stiffnessMatrix, loadVector);

            std::vector<double> nodes = FillNodes(minimum, matrixSize, step / EElementTypeToInt[options.Type]);
            std::vector<double> displacementsReal = FillRealValues(nodes, RealSolve);
            std::vector<double> errors = FillErrors(displacements.GetColumn(0), displacementsReal);

            double maxError = (*std::max_element(errors.begin(), errors.end()));

            std::cout << ke_type << " " << size << std::endl;
            std::cout << "Maximal error: " << maxError << std::endl;

            {
                std::string filename = ke_type + "_" + std::to_string(size);
                TPlotter<double> graphics(filename);
                graphics.SetXValues(nodes);
                graphics.EmplaceGraphic<double>("displacements", displacements.GetColumn(0));
                graphics.Plot();
            }

            {
                double d2 = (displacements[2][0] - displacements[1][0]) / step;
                double d1 = (displacements[1][0] - displacements[0][0]) / step;
                double value = d1 - (d2 - d1) / 2;
                std::cout << "Absolute error: " << std::fabs(value - 10) << std::endl;
                std::cout << "Relative error: " << (std::fabs(value - 10) / std::fabs(value)) * 100 << "%" << std::endl;
            }

            {
                double d2 = (displacements[matrixSize - 1][0] - displacements[matrixSize - 2][0]) / step;
                double d1 = (displacements[matrixSize - 2][0] - displacements[matrixSize - 3][0]) / step;
                double value = d1 + (d2 - d1) / 2;
                std::cout << "Absolute error: " << std::fabs(value - displacements[matrixSize - 1][0]) << std::endl;
                std::cout << "Relative error: " << (std::fabs(value - displacements[matrixSize - 1][0]) / std::fabs(value)) * 100 << "%" << std::endl;
            }

            MoveFilesToFolder(".png", "images");
            MoveFilesToFolder(".csv", "results");
        }
    }
    return 0;
}
