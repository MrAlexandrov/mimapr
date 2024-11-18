#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include "matrix.hpp"

using namespace NMatrix;

enum element_type {unknown, linear, cubic = 3};
enum restrict_grade {first, second, third};


// Структура для передачи опций программе
typedef struct {
    element_type type;
    int elemAmount;
    int help = 0;
} opt;

// Структура для описания граничных условий
typedef struct {
    restrict_grade grade;
    double pos;
    double val;
} restriction;


void printBold(const std::string& text) {
    std::cout << "\x1B[1m" << text << "\x1B[0m";
}

void printUsage() {
    using std::cout;
    using std::endl;

    // Заголовок раздела
    cout << "Usage:\n"
         << "./mke [-l | -c] [-s <SIZE>]\n\n";

    // Описание опций
    cout << "Options:\n";
    
    printBold("\t-h, --help");
    cout << "\n\t\tdisplay help.\n";

    printBold("\t-s, --size");
    cout << "\n\t\tset amount of elements (default is 20).\n";

    printBold("\t-l, --linear");
    cout << "\n\t\tUse linear equations for elements (used by default).\n";

    printBold("\t-c, --cubic");
    cout << "\n\t\tUse cubic equations for elements.\n";

    cout << endl;
}

// Обработка опций программы
opt getOptions(const int argc, char *argv[]) {
    opt options;
    options.type = unknown;
    options.elemAmount = 0;

    const struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"linear", no_argument, NULL, 'l'},
        {"cubic", no_argument, NULL, 'c'},
        {"size", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};
    int long_optind;
    char c;
    while ((c = getopt_long(argc, argv, "hlcs:", long_options, &long_optind)) != -1) {
        switch (c) {
        case 'h':
            printUsage();
            options.help = 1;
            return options;
        case 'l':
            options.type = linear;
            break;
        case 'c':
            options.type = cubic;
            break;
        case 's':
            options.elemAmount = atoi(optarg);
            break;
        default:
            break;
        }
    }

    if (!options.elemAmount) {
        std::cout << "Element count was set to default(20)." << std::endl;
        options.elemAmount = 20;
    }

    if (options.type == unknown) {
        std::cout << "Function type was set to default(linear)." << std::endl;
        options.type = linear;
    }

    return options;
}


// Заполнение СЛАУ линейными КЭ
int makeLinearSLAU(TMatrix<> &resultMatrix, 
                   TMatrix<> &resultVector,
                   double L, 
                   double a, double b, double c, double d, 
                   int i) 
{
    resultMatrix[i][i]              += -a / L - b / 2. + c * L / 3.;
    resultMatrix[i][i + 1]          +=  a / L + b / 2. + c * L / 6.;
    resultMatrix[i + 1][i]          +=  a / L - b / 2. + c * L / 6.;
    resultMatrix[i + 1][i + 1]      += -a / L + b / 2. + c * L / 3.;

    resultVector[i][0]              += -d * L / 2.;
    resultVector[i + 1][0]          += -d * L / 2.;
    return 0;
}

// Заполнение СЛАУ кубическими КЭ
int makeCubicSLAU(TMatrix<> &resultMatrix, 
                  TMatrix<> &resultVector,
                  double L, 
                  double a, double b, double c, double d, 
                  int i) 
{
    resultMatrix[i][i]              += -a *  37 / (10 * L) - b / 2       + c * 8  * L / 105;
    resultMatrix[i][i + 1]          +=  a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[i][i + 2]          += -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[i][i + 3]          +=  a *  13 / (40 * L) + b * 7  / 80 + c * 19 * L / 1680;

    resultMatrix[i + 1][i]          +=  a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[i + 1][i + 1]      += -a *  54 / (5  * L)                + c * 27 * L / 70;
    resultMatrix[i + 1][i + 2]      +=  a * 297 / (40 * L) + b * 81 / 80 - c * 27 * L / 560;
    resultMatrix[i + 1][i + 3]      += -a *  27 / (20 * L) - b * 3  / 10 - c * 3  * L / 140;

    resultMatrix[i + 2][i]          += -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[i + 2][i + 1]      +=  a * 297 / (40 * L) - b * 81 / 80 - c * 27 * L / 560;
    resultMatrix[i + 2][i + 2]      += -a *  54 / (5  * L)                + c * 27 * L / 70;
    resultMatrix[i + 2][i + 3]      +=  a * 189 / (40 * L) + b * 57 / 80 + c * 33 * L / 560;

    resultMatrix[i + 3][i]          +=  a *  13 / (40 * L) - b * 7  / 80 + c * 19 * L / 1680;
    resultMatrix[i + 3][i + 1]      += -a *  27 / (20 * L) + b * 3  / 10 - c * 3  * L / 140;
    resultMatrix[i + 3][i + 2]      +=  a * 189 / (40 * L) - b * 57 / 80 + c * 33 * L / 560;
    resultMatrix[i + 3][i + 3]      += -a *  37 / (10 * L) + b / 2       + c * 8  * L / 105;

    resultVector[i][0]              += -d * L / 8;
    resultVector[i + 1][0]          += -d * 3 * L / 8;
    resultVector[i + 2][0]          += -d * 3 * L / 8;
    resultVector[i + 3][0]          += -d * L / 8;
    return 0;
}


// Применение граничных условий
int restrictType1(TMatrix<>& matrix, TMatrix<>& vector, int a, int row, double value) {
    for (int col = 0, end = matrix.cols(); col < end; ++col) {
        matrix[row][col] = 0;
    }
    matrix[row][row] = 1;
    vector[row][0] = value;
    return 0;
}

int restrictType2(TMatrix<>& matrix, TMatrix<>& vector, int a, int row, double value) {
    double coef = 0;
    if (row == 0) {
        coef = 1;
    } else if (row == vector.rows() - 1) {
        coef = -1;
    } else {
        return 0;
    }
    vector[row][0] += a * coef * value;
    return 0;
}

int restrictType3(TMatrix<>& matrix, TMatrix<>& vector, int a, int row, double value) {
    if (row == 0) {
        matrix[row][row] += value * -a;
    } else if (row == matrix.rows() - 1) {
        matrix[row][row] += value * a;
    }
    return 0;
}


int (*restrict[3])(TMatrix<>& matrix, TMatrix<>& vector, int a, int line, double value) = {
    restrictType1, restrictType2, restrictType3};


int plot(const std::string& filename, double left, double right) {
    // Открываем канал для общения с gnuplot
    FILE* gp = popen("gnuplot -persist", "w");
    if (!gp) {
        std::cerr << "Gnuplot не найден\n";
        return 1;
    }

    // Формируем команды для gnuplot
    std::ostringstream commands;
    commands << "set term png size 640, 480\n";
    commands << "set output '" << filename << ".png'\n";
    commands << "set title 'Displacements'\n";
    commands << "set xlabel 'X'\n";
    commands << "set ylabel 'U(X)'\n";
    commands << "set xrange [" << left << ":" << right << "]\n";
    commands << "plot '" << filename << ".txt' using 1:3 with lines title 'MKE', ";
    commands << "'" << filename << ".txt' using 1:2 with lines title 'Real'\n";

    // Записываем команды в поток gnuplot
    std::string cmd_str = commands.str();
    fwrite(cmd_str.c_str(), sizeof(char), cmd_str.size(), gp);

    // Закрываем канал
    pclose(gp);
    return 0;
}

// Решение СЛАУ методом Гаусса
int solveSLAU(TMatrix<> &retXVector, TMatrix<> inAMatrix, TMatrix<> inBVector, int size) {
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            double d = inAMatrix[j][i] / inAMatrix[i][i];
            for (int k = i; k < size; k++)
                inAMatrix[j][k] -= d * inAMatrix[i][k];
            inBVector[j][0] -= d * inBVector[i][0];
        }
    }

    for (int i = size - 1; i >= 0; --i) {
        double rSum = 0;
        for (int j = i + 1; j < size; ++j) {
            rSum += retXVector[j][0] * inAMatrix[i][j];
        }
        retXVector[i][0] = (inBVector[i][0] - rSum) / inAMatrix[i][i];
    }
    return 0;
}

// Точное решение
double realSolve(double x) {
	double C1 = -(5 * exp(8) * (-3 + 2 * exp(6))) / (1 + exp(12));
	double C2 = 5 * (2 + 3 * exp(6)) / (exp(2) * (1 + exp(12)));
	return exp(-x) * C1 + exp(x) * C2 - 10;
}


int main(int argc, char **argv) {
    // Обработка опций программы
    const opt settings = getOptions(argc, argv);
    if (settings.help) {
        return 0;
    }

    // Данные по варианту
    const double a = 1, b = 0, c = -1, d = -10;
    restriction lower = {second, 2, 10},
                upper = {first, 8, 5};

    // Вычисление необходимых постоянных
    const int size = settings.elemAmount * settings.type + 1;
    const double step = (upper.pos - lower.pos) / (settings.elemAmount);

    // Используемые структуры данных
    std::fstream out;
    TMatrix<> stiffnessMatrix(size, size, 0);
    TMatrix<> loadVector(size, 1);
    TMatrix<> displacements(size, 1);
    std::vector<double> errors;
    std::vector<double> nodes;
    std::vector<double> displacementsReal;

    // Заполняем матрицу жесткости и вектор нагрузок
    if (settings.type == linear) {
        for (int i = 0; i < size - 1; ++i)
            makeLinearSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
    } else if (settings.type == cubic) {
        for (int i=0; i < size - 3; i+=3)
            makeCubicSLAU(stiffnessMatrix, loadVector, step, a, b, c, d, i);
    } else {
        std::cout << "Some error occured, equation type not specified." << std::endl;
        return -1;
    }
    
    // Применяем граничные условия
    restrict[lower.grade](stiffnessMatrix, loadVector, a, 0, lower.val);
    restrict[upper.grade](stiffnessMatrix, loadVector, a, size - 1, upper.val);

    // Решаем СЛАУ методом Гаусса
    solveSLAU(displacements, stiffnessMatrix, loadVector, size);

    // Вычисляем погрешности и заполняем массивы для вывода графиков
    double node = lower.pos;
    for (int i = 0; i < size; i++) {
        double realU = realSolve(node);
        double error = fabs(displacements[i][0] - realU);

        errors.push_back(error);
        nodes.push_back(node);
        displacementsReal.push_back(realU);

        std::cout << "u(" << node << ") = " << displacements[i][0] << std::endl;
        node += (step) / settings.type;
    }
    double maxError = *std::max_element(errors.begin(), errors.end());

    // Выводим наибольшую погрешность
    std::cout << std::endl << "Max error: " << maxError << "\n\n";

    // Создаём файл с результатами рассчетов
    std::string filename = {"lc"[settings.type / 2], '_'};
    filename += std::to_string((size - 1) / settings.type);
    out.open(filename + ".txt", std::ios::out);
    if (out.is_open()) {
        for (int i = 0; i < size; i++)
            out << nodes[i] << '\t' << displacementsReal[i] << '\t' << displacements[i][0] << '\t' << errors[i] << std::endl;
        // Выводим графики через Gnuplot
        plot(filename.c_str(), lower.pos, upper.pos);
    } else
        std::cerr << "Output file didn't open." << std::endl;
    return 0;
}
