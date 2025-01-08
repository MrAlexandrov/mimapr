#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>

#include "initialize.hpp"
#include "matrix.hpp"
#include "options.hpp"
#include "plotter.hpp"

namespace {

using namespace NMatrix;
using namespace NPlotter;
using namespace NOptions;

void MakeLinearSLAU(TMatrix<>&          resultMatrix, 
                    TMatrix<>&          resultVector, 
                    const long double   L,
                    const int           index)
{
    resultMatrix[index][index]          += -A / L   - B / 2 + C * L / 3;
    resultMatrix[index][index + 1]      +=  A / L   + B / 2 + C * L / 6;
    resultMatrix[index + 1][index]      +=  A / L   - B / 2 + C * L / 6;
    resultMatrix[index + 1][index + 1]  += -A / L   + B / 2 + C * L / 3;

    resultVector[index][0]              += -D * L   / 2;
    resultVector[index + 1][0]          += -D * L   / 2;
}

void MakeCubicSLAU(TMatrix<>&           resultMatrix, 
                   TMatrix<>&           resultVector, 
                   const long double    L,
                   const int            index) 
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

void InitializeSLAU(TMatrix<>&              stiffnessMatrix,
                    TMatrix<>&              loadVector,
                    const EElementType&     type,
                    const int               size,
                    const long double       step) 
{
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

void ApplyRestriction(TMatrix<>&            matrix, 
                      TMatrix<>&            vector, 
                      const TRestriction&   restriction, 
                      const int             row) 
{
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

void SolveSLAU(TMatrix<>& xVector, TMatrix<>& aMatrix, TMatrix<>& bVector) {
    int size = aMatrix.Rows();
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

inline long double RealSolve(const long double x) {
	long double C1 = -(5 * expl(8) * (-3 + 2 * expl(6))) / (1 + expl(12));
	long double C2 = 5 * (2 + 3 * expl(6)) / (expl(2) * (1 + expl(12)));
	return expl(-x) * C1 + expl(x) * C2 - 10;
}

inline long double RealSolveAdditional(const long double x) {
    return static_cast<long double>(5) * (expl(8 - x) + 2 * expl(x - 2) + expl(x + 4) - 2);
}

std::vector<long double> FillNodes(
                                    long double minimum,
                                    long double size, 
                                    long double add) 
{
    std::vector<long double> result;
    long double current = minimum;
    for (int i = 0; i < size; ++i, current += add) {
        result.push_back(current);
    }
    return result;
}

std::vector<long double> FillRealValues(const std::vector<long double>& positions, const auto& function) {
    std::vector<long double> result;
    for (const auto& x : positions) {
        result.emplace_back(function(x));
    }
    return result;
}

std::vector<long double> FillErrors(const std::vector<long double>& solution,
                       const std::vector<long double>& realValues) 
{
    assert(solution.size() == realValues.size());
    std::vector<long double> result;
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

} // namespace

int main(int argc, char* argv[]) {
    const std::optional<TOptions> options = GetOptions(argc, argv);
    if (!options.has_value() || options->Help) {
        return 0;
    }

    const int size = options->ElementsAmount * EElementTypeToInt[options->Type] + 1;
    const long double step = static_cast<long double>(maximum - minimum) / options->ElementsAmount;

    TMatrix<> stiffnessMatrix(size, size, 0.0);
    TMatrix<> loadVector(size, 1, 0.0);
    TMatrix<> displacements(size, 1, 0.0);

    InitializeSLAU(stiffnessMatrix, loadVector, options->Type, size, step);

    auto getRowByPosition = [&](long double position) -> int {
        assert(minimum <= position && position <= maximum);
        long double shift = position - minimum;
        shift *= EElementTypeToInt[options->Type];
        int row = std::round(shift / step);
        row = std::max(0, row);
        row = std::min(size - 1, row);
        return row;
    };

    for (const auto& current : restrictions) {
        ApplyRestriction(
            stiffnessMatrix,
            loadVector,
            current, 
            getRowByPosition(current.Position)
        );
    }

    SolveSLAU(displacements, stiffnessMatrix, loadVector);

    long double add = step / EElementTypeToInt[options->Type];
    std::vector<long double> nodes = FillNodes(minimum, size, add);
    std::vector<long double> displacementsReal = FillRealValues(nodes, RealSolve);
    std::vector<long double> errors = FillErrors(displacements.GetColumn(0), displacementsReal);

    // std::cout << nodes << std::endl;
    // std::cout << displacementsReal << std::endl;
    // std::cout << displacements << std::endl;
    // std::cout << errors << std::endl;
    // return 0;

    long double maxError = (*std::max_element(errors.begin(), errors.end()));

    std::cout << EElementTypeToString[options->Type] << ' ';
    std::cout << options->ElementsAmount << std::endl;
    std::cout << "Maximal error: " << maxError << std::endl;

    {
        const std::string filename = 
            EElementTypeToString[options->Type]
            + std::string("_") 
            + std::to_string(options->ElementsAmount);

        TPlotter graphics(filename);
        graphics.SetXRangeLeft(minimum);
        graphics.SetXRangeRight(maximum);
        graphics.SetXValues(nodes);
        graphics.AddGraphic("displacements", displacements.GetColumn(0));
        // graphics.AddGraphic("displacements real", displacementsReal);
        graphics.Plot();

        {
        long double d2 = (displacements[2][0] - displacements[1][0]) / add;
        long double d1 = (displacements[1][0] - displacements[0][0]) / add;
        // std::cout << d2 << std::endl;
        // std::cout << d1 << std::endl;
        long double value = d1 - (d2 - d1) / 2;
        // std::cout << "Approximate du/dx(2): " << value << std::endl;
        std::cout << "Absolute error: " << std::fabs(value - 10) << std::endl;
        std::cout << "Relative error: " << (std::fabs(value - 10) / std::fabs(value)) * 100 << "%" << std::endl;
        }

        {
        long double d2 = (displacements[size - 1][0] - displacements[size - 2][0]) / add;
        long double d1 = (displacements[size - 2][0] - displacements[size - 3][0]) / add;
        // std::cout << d2 << std::endl;
        // std::cout << d1 << std::endl;
        long double value = d1 + (d2 - d1) / 2;
        // std::cout << "Approximate du/dx(8)=u, U(8): " << value << ", " << displacements[size - 1][0] << std::endl;
        std::cout << "Absolute error: " << std::fabs(value - displacements[size - 1][0]) << std::endl;
        std::cout << "Relative error: " << (std::fabs(value - displacements[size - 1][0]) / std::fabs(value)) * 100 << "%" << std::endl;
        }
    }

    const int derivativeSize = size - 1;
    std::vector<long double> derivativeNodes;
    std::vector<long double> derivativeDisplacements;

    for (int i = 0; i < derivativeSize; ++i) {
        derivativeNodes.emplace_back(nodes[i] + (nodes[i + 1] - nodes[i]) / 2);
        derivativeDisplacements.emplace_back((displacements[i + 1][0] - displacements[i][0]) / step);
    }

    {
        const std::string filename =
            std::string("DERIVATIVE_")
            + EElementTypeToString[options->Type]
            + std::string("_") 
            + std::to_string(options->ElementsAmount);

        TPlotter derivative(filename);
        derivative.SetXRangeLeft(minimum);
        derivative.SetXRangeRight(maximum);
        derivative.SetXValues(derivativeNodes);
        derivative.AddGraphic("derivativeDisplacements", derivativeDisplacements);
        derivative.Plot();
    }
    return 0;
}
