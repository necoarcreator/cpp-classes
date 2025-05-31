import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def find_min_max(folder_num, folder_ana, fict):
    """
    Находит глобальные минимумы и максимумы для каждой величины (p, rho, v, e).
    """
    min_max = {"p": [float('inf'), float('-inf')],
               "rho": [float('inf'), float('-inf')],
               "v": [float('inf'), float('-inf')],
               "e": [float('inf'), float('-inf')]}

    def update_min_max(values, key):
        min_max[key][0] = min(min_max[key][0], np.min(values[fict:-fict]))
        min_max[key][1] = max(min_max[key][1], np.max(values[fict:-fict]))

    num_files = sorted([f for f in os.listdir(folder_num) if f.endswith('.csv')])
    ana_files = sorted([f for f in os.listdir(folder_ana) if f.endswith('.csv')])

    for num_file, ana_file in zip(num_files, ana_files):
        # Численное решение
        df_num = pd.read_csv(os.path.join(folder_num, num_file), sep=",", header=1, skiprows=0,
                             names=["x", "v", "rho", "e", "p"])
        update_min_max(pd.to_numeric(df_num["p"]), "p")
        update_min_max(pd.to_numeric(df_num["rho"]), "rho")
        update_min_max(pd.to_numeric(df_num["v"]), "v")
        update_min_max(pd.to_numeric(df_num["e"]), "e")

        # Аналитическое решение
        df_ana = pd.read_csv(os.path.join(folder_ana, ana_file), sep=",", header=1, skiprows=0,
                             names=["x", "v", "rho", "e", "p"])
        update_min_max(pd.to_numeric(df_ana["p"]), "p")
        update_min_max(pd.to_numeric(df_ana["rho"]), "rho")
        update_min_max(pd.to_numeric(df_ana["v"]), "v")
        update_min_max(pd.to_numeric(df_ana["e"]), "e")

    return min_max

def Plot(x_num, y_num, x_ana, y_ana, time, output_path, file_index, min_max):
    """
    Построение графиков численного и аналитического решений с одинаковыми масштабами.
    """
    name = [[f'p', f'rho'], [f'v', f'e']]
    keys = [["p", "rho"], ["v", "e"]]
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

    for i in range(2):
        for j in range(2):
            key = keys[i][j]
            # Численное решение
            axes[i][j].scatter(x=x_num[fict:Nx - fict], y=y_num[i][j][fict:Nx - fict], marker='o', c='b', label='Численное')
            # Аналитическое решение
            axes[i][j].plot(x_ana, y_ana[i][j], linestyle='--', c='r', label='Аналитическое')

            axes[i][j].set_title(f'${name[i][j]}$, time = {time}')
            axes[i][j].set_xlabel('$x$')
            axes[i][j].set_ylabel(f'${name[i][j]}$')
            axes[i][j].set_xlim([0, 1.0])
            axes[i][j].set_ylim(min_max[key][0], min_max[key][1])
            axes[i][j].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, f'plot_{file_index:03d}.png'))
    plt.close()

# Запрос папок
num_folder = input("Введите номер папки с численным решением (1, 2, 3...): ").strip()
ana_folder = input("Введите номер папки с аналитическим решением (1, 2, 3...): ").strip()

# Формирование путей
num_folder = f"results{num_folder}"
ana_folder = f"results{ana_folder}"

# Проверка существования папок
if not os.path.isdir(num_folder) or not os.path.isdir(ana_folder):
    print("Одна из указанных папок не существует.")
    exit(1)

# Получение списка файлов .csv
num_files = sorted([f for f in os.listdir(num_folder) if f.endswith('.csv')])
ana_files = sorted([f for f in os.listdir(ana_folder) if f.endswith('.csv')])

if not num_files or not ana_files:
    print("В одной из папок нет файлов с расширением .csv.")
    exit(1)

# Проверка соответствия количества файлов
if len(num_files) != len(ana_files):
    print("Число файлов в папках не совпадает.")
    exit(1)

# Найти минимумы и максимумы для всех величин
fict = 1  # Параметр для обрезки данных
min_max = find_min_max(num_folder, ana_folder, fict)

# Построение графиков
for index, (num_file, ana_file) in enumerate(zip(num_files, ana_files)):
    # Численное решение
    with open(os.path.join(num_folder, num_file), 'r') as file:
        time = float(file.readline().strip())
    df_num = pd.read_csv(os.path.join(num_folder, num_file), sep=",", header=1, skiprows=0,
                         names=["x", "v", "rho", "e", "p"])
    x_num = pd.to_numeric(df_num["x"]).to_numpy()
    v_num = pd.to_numeric(df_num["v"]).to_numpy()
    rho_num = pd.to_numeric(df_num["rho"]).to_numpy()
    e_num = pd.to_numeric(df_num["e"]).to_numpy()
    p_num = pd.to_numeric(df_num["p"]).to_numpy()
    y_num = [[p_num, rho_num], [v_num, e_num]]

    # Аналитическое решение
    with open(os.path.join(ana_folder, ana_file), 'r') as file:
        time = float(file.readline().strip())
    df_ana = pd.read_csv(os.path.join(ana_folder, ana_file), sep=",", header=1, skiprows=0,
                         names=["x", "v", "rho", "e", "p"])
    x_ana = pd.to_numeric(df_ana["x"]).to_numpy()
    v_ana = pd.to_numeric(df_ana["v"]).to_numpy()
    rho_ana = pd.to_numeric(df_ana["rho"]).to_numpy()
    e_ana = pd.to_numeric(df_ana["e"]).to_numpy()
    p_ana = pd.to_numeric(df_ana["p"]).to_numpy()
    y_ana = [[p_ana, rho_ana], [v_ana, e_ana]]

    Nx = np.size(x_num)

    # Построение графика
    Plot(x_num, y_num, x_ana, y_ana, time, num_folder, index + 1, min_max)

print(f"Графики сохранены в папке: {num_folder}")
print("Для создания анимации используйте команду:")
print(f'ffmpeg -r 5 -i "{num_folder}/plot_%03d.png" animation.mp4')
