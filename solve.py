import subprocess
import os


def run_prog(code_folder):
    build_result = subprocess.run(['make', 'build'], cwd=code_folder)
    if build_result.returncode != 0:
        raise RuntimeError("Сборка завершилась с ошибкой.")

    results_table = {}
    sizes = [1, 3, 20, 40, 100, 200, 400, 800, 1600, 3200]
    ke_types = ['linear', 'cubic']

    program_path = os.path.join(code_folder, 'build', 'main')

    for size in sizes:
        row = {}
        for ke_type in ke_types:
            if (ke_type == 'linear' and size == 1) or (ke_type == 'cubic' and size == 3):
                row[ke_type] = '-'
                continue

            result = subprocess.run(
                [program_path, f'--{ke_type}', '-s', str(size)],
                capture_output=True, text=True
            )

            if result.stdout:
                print(result.stdout)

            if result.returncode != 0:
                raise RuntimeError(f"Ошибка выполнения {program_path} с параметрами -{ke_type} -s {size}")

    img_folder = os.path.join(code_folder, 'images')
    res_folder = os.path.join(code_folder, 'results')

    os.makedirs(img_folder, exist_ok=True)
    os.makedirs(res_folder, exist_ok=True)

    for file in os.listdir(code_folder):
        if file.endswith('.png'):
            os.rename(os.path.join(code_folder, file), os.path.join(img_folder, file))

    for file in os.listdir(code_folder):
        if file.endswith('.csv'):
            os.rename(os.path.join(code_folder, file), os.path.join(res_folder, file))

    return


if __name__ == "__main__":
    try:
        project_path = "."
        results = run_prog(project_path)
    except RuntimeError as error:
        print(f"Ошибка: {error}")