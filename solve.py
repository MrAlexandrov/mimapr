import subprocess
import os


def run_prog(code_folder):
    # Шаг 1: Сборка программы с использованием make
    build_result = subprocess.run(['make', 'build'], cwd=code_folder)
    if build_result.returncode != 0:
        raise RuntimeError("Сборка завершилась с ошибкой.")

    # Шаг 2: Создание таблицы для хранения результатов
    results_table = {}
    sizes = [1, 3, 20, 40]
    ke_types = ['l', 'c']  # Типы данных или вычислений: 'l' и 'c'

    program_path = os.path.join(code_folder, 'build', 'main')

    for size in sizes:
        row = {}
        for ke_type in ke_types:
            # Пропуск комбинаций 'l' с size=1 и 'c' с size=3
            if (ke_type == 'l' and size == 1) or (ke_type == 'c' and size == 3):
                row[ke_type] = '-'  # Указать, что этот случай пропускается
                continue

            # Запуск программы ./build/main с текущими параметрами
            result = subprocess.run(
                [program_path, f'-{ke_type}', '-s', str(size)],
                capture_output=True, text=True
            )

            # Проверка на успешное выполнение
            if result.returncode != 0:
                raise RuntimeError(f"Ошибка выполнения {program_path} с параметрами -{ke_type} -s {size}")

    # Шаг 3: Организация выходных файлов
    img_folder = os.path.join(code_folder, 'images')
    res_folder = os.path.join(code_folder, 'results')

    # Создание папок для изображений и результатов, если они ещё не существуют
    os.makedirs(img_folder, exist_ok=True)
    os.makedirs(res_folder, exist_ok=True)

    # Перемещение всех изображений (*.png) в папку images
    for file in os.listdir(code_folder):
        if file.endswith('.png'):
            os.rename(os.path.join(code_folder, file), os.path.join(img_folder, file))

    # Перемещение всех текстовых файлов (*.txt) в папку results
    for file in os.listdir(code_folder):
        if file.endswith('.txt') and file != "CMakeLists.txt":
            os.rename(os.path.join(code_folder, file), os.path.join(res_folder, file))

    return


if __name__ == "__main__":
    try:
        project_path = "."  # Корневая директория проекта
        results = run_prog(project_path)
    except RuntimeError as e:
        print(f"Ошибка: {e}")