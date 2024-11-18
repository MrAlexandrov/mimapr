#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <optional>
#include "matrix.hpp"

using namespace NMatrix;

// Тип элемента
enum class ElementType { Unknown, Linear, Cubic };

// Тип ограничения
enum class RestrictGrade { First, Second, Third };

// Настройки программы
struct Options {
    ElementType type = ElementType::Unknown;
    int elemAmount = 20;  // По умолчанию
    bool help = false;
};

// Граничное условие
struct Restriction {
    RestrictGrade grade;
    double pos;
    double val;
};

// Вывод жирного текста в консоль
void printBold(const std::string& text) {
    std::cout << "\x1B[1m" << text << "\x1B[0m";
}

// Вывод инструкции по использованию программы
void printUsage() {
    std::cout << "Usage:\n"
              << "./mke [-l | -c] [-s <SIZE>]\n\n"
              << "Options:\n";

    printBold("\t-h, --help");
    std::cout << "\n\t\tdisplay help.\n";

    printBold("\t-s, --size");
    std::cout << "\n\t\tset amount of elements (default is 20).\n";

    printBold("\t-l, --linear");
    std::cout << "\n\t\tUse linear equations for elements (used by default).\n";

    printBold("\t-c, --cubic");
    std::cout << "\n\t\tUse cubic equations for elements.\n";
}

// Получение настроек программы
std::optional<Options> getOptions(int argc, char* argv[]) {
    Options options;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage();
            options.help = true;
            return options;
        } else if (arg == "-l" || arg == "--linear") {
            options.type = ElementType::Linear;
        } else if (arg == "-c" || arg == "--cubic") {
            options.type = ElementType::Cubic;
        } else if (arg == "-s" || arg == "--size") {
            if (i + 1 < argc) {
                options.elemAmount = std::stoi(argv[++i]);
            } else {
                std::cerr << "Missing value for option -s.\n";
                return std::nullopt;
            }
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            return std::nullopt;
        }
    }

    if (options.type == ElementType::Unknown) {
        std::cout << "Element type not specified. Using default: Linear.\n";
        options.type = ElementType::Linear;
    }

    return options;
}

// Создание матрицы жесткости (линейные элементы)
void makeLinearSLAU(TMatrix<>& resultMatrix, TMatrix<>& resultVector, double L, double a, double b, double c, double d, int i) {
    resultMatrix[i][i]          += -a / L   - b / 2.0 + c * L / 3.0;
    resultMatrix[i][i + 1]      += a / L    + b / 2.0 + c * L / 6.0;
    resultMatrix[i + 1][i]      += a / L    - b / 2.0 + c * L / 6.0;
    resultMatrix[i + 1][i + 1]  += -a / L   + b / 2.0 + c * L / 3.0;

    resultVector[i][0]          += -d * L   / 2.0;
    resultVector[i + 1][0]      += -d * L   / 2.0;
}

// Создание матрицы жесткости (кубические элементы)
void makeCubicSLAU(TMatrix<>& resultMatrix, 
                   TMatrix<>& resultVector, 
                   double L, 
                   double a, 
                   double b, 
                   double c, 
                   double d, 
                   int i) {
    // Матрица жёсткости (stiffness matrix)
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

    // Вектор нагрузки (load vector)
    resultVector[i][0]          += -d * L / 8;
    resultVector[i + 1][0]      += -d * 3 * L / 8;
    resultVector[i + 2][0]      += -d * 3 * L / 8;
    resultVector[i + 3][0]      += -d * L / 8;
}

// Применение граничных условий
void applyRestrictions(TMatrix<>& matrix, TMatrix<>& vector, Restriction restriction, double coef) {
    int row = (restriction.grade == RestrictGrade::First ? 0 : matrix.rows() - 1);
    if (restriction.grade == RestrictGrade::Second) {
        vector[row][0] += coef * restriction.val;
    } else {
        matrix[row][row] = 1.0;
        vector[row][0] = restriction.val;
    }
}

// Решение СЛАУ методом Гаусса
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

// Главная функция программы
int main(int argc, char* argv[]) {
    auto opt = getOptions(argc, argv);
    if (!opt.has_value() || opt->help) {
        return 0;
    }

    const double a = 1, b = 0, c = -1, d = -10;
    Restriction lower = {RestrictGrade::Second, 2, 10};
    Restriction upper = {RestrictGrade::First, 8, 5};

    int size = opt->elemAmount * (opt->type == ElementType::Linear ? 1 : 3) + 1;
    double step = (upper.pos - lower.pos) / opt->elemAmount;

    TMatrix<> stiffnessMatrix(size, size, 0.0);
    TMatrix<> loadVector(size, 1, 0.0);
    TMatrix<> displacements(size, 1, 0.0);

    // Заполнение матриц
    for (int i = 0; i < size - 1; ++i) {
        makeLinearSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
    }

    applyRestrictions(stiffnessMatrix, loadVector, lower, a);
    applyRestrictions(stiffnessMatrix, loadVector, upper, a);

    solveSLAU(displacements, stiffnessMatrix, loadVector);

    for (int i = 0; i < size; ++i) {
        std::cout << "U(" << i << ") = " << displacements[i][0] << "\n";
    }

    return 0;
}
