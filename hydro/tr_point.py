from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob
#import mplcyberpunk
#plt.style.use("cyberpunk")

# Настройки визуализации
PLOT_SETTINGS = {
    'cell': {
        'edgecolor': 'k',
        'linewidth': 0.3,
        'antialiased': False,
        'shading': 'nearest'
    },
    'grid': {
        'color': 'black',
        'linestyle': '--',
        'linewidth': 0.25,
        'alpha': 0
    },
    'figure': {
        'size': (20, 12),
        'dpi': 100
    }
}

def read_data(file_path):
    """Read and parse data file into grid structure"""
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        raw_data = [line.strip().rstrip(';').split(',') for line in f if line.strip()]
    
    data = np.array([[float(x) for x in row] for row in raw_data])
    
    # Создание регулярной сетки
    x_vals = np.unique(data[:, 0])
    y_vals = np.unique(data[:, 1])
    X, Y = np.meshgrid(x_vals, y_vals)
    
    # Инициализация параметров
    params = {
        'X': X,
        'Y': Y,
        'p': np.full_like(X, np.nan),
        'vx': np.full_like(X, np.nan),
        'vy': np.full_like(X, np.nan),
        'r': np.full_like(X, np.nan),
        'e': np.full_like(X, np.nan)
    }
    
    # Заполнение данных
    coord_map = {(x, y): idx for idx, (x, y) in enumerate(zip(data[:,0], data[:,1]))}
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            key = (X[i,j], Y[i,j])
            if key in coord_map:
                idx = coord_map[key]
                params['p'][i,j] = data[idx,2]
                params['vx'][i,j] = data[idx,3]
                params['vy'][i,j] = data[idx,4]
                params['r'][i,j] = data[idx,5]
                params['e'][i,j] = data[idx,6]
    
    # Статистика для файла
    stats = {}
    for param in ['p', 'vx', 'vy', 'r', 'e']:
        valid_vals = params[param][~np.isnan(params[param])]
        stats[param] = (np.min(valid_vals), np.max(valid_vals)) if valid_vals.size > 0 else (np.inf, -np.inf)
    
    return params, stats

def get_global_stats(comm, local_stats):
    """Вычисление глобальных минимумов и максимумов"""
    global_stats = {}
    for param in ['p', 'vx', 'vy', 'r', 'e']:
        # Сбор локальных значений
        local_mins = comm.gather(local_stats[param][0], root=0)
        local_maxs = comm.gather(local_stats[param][1], root=0)
        
        # Вычисление глобальных значений
        if comm.Get_rank() == 0:
            global_min = np.min([m for m in local_mins if m != np.inf])
            global_max = np.max([m for m in local_maxs if m != -np.inf])
        else:
            global_min = global_max = None
        
        # Рассылка всем процессам
        global_min = comm.bcast(global_min, root=0)
        global_max = comm.bcast(global_max, root=0)
        global_stats[param] = (global_min, global_max)
    
    return global_stats

def create_plots(data, output_path, global_stats):
    """Генерация графиков с фиксированными шкалами"""
    plt.figure(figsize=PLOT_SETTINGS['figure']['size'], 
              dpi=PLOT_SETTINGS['figure']['dpi'])
    
    # Общие настройки
    plt.gca().set_aspect('equal')
    
    # Создание графиков
    def create_subplot(pos, data, param, label, title):
        plt.subplot(2, 3, pos)

        mesh = plt.pcolormesh(data['X'], data['Y'], data[param],
                            #vmin=global_stats[param][0],
                            #vmax=global_stats[param][1],
                            cmap='viridis',
                            **PLOT_SETTINGS['cell'])
        plt.colorbar(mesh, label=label)
        plt.title(title)
        #plt.grid(**PLOT_SETTINGS['grid'])
        plt.xlabel('X')
        plt.ylabel('Y')
    
    # Скорость по X
    create_subplot(1, data, 'vx', 'Vx', 'X-Velocity')
    
    # Скорость по Y
    create_subplot(2, data, 'vy', 'Vy', 'Y-Velocity')
    
    # Модуль скорости
    plt.subplot(2, 3, 3)
    speed = np.sqrt(data['vx']**2 + data['vy']**2)
    #speed_min = np.sqrt(global_stats['vx'][0]**2 + global_stats['vy'][0]**2)
    #speed_max = np.sqrt(global_stats['vx'][1]**2 + global_stats['vy'][1]**2)
    mesh = plt.pcolormesh(data['X'], data['Y'], speed,
                        #vmin=speed_min, vmax=speed_max,
                        cmap='viridis',
                        **PLOT_SETTINGS['cell'])
    plt.colorbar(mesh, label='Speed')
    plt.title('Speed Magnitude')
    #plt.grid(**PLOT_SETTINGS['grid'])
    plt.xlabel('X')
    plt.ylabel('Y')
    
    # Давление
    create_subplot(4, data, 'p', 'Pressure', 'Pressure')
    
    # Плотность
    create_subplot(5, data, 'r', 'Density', 'Density')
    
    # Энергия
    create_subplot(6, data, 'e', 'Energy', 'Energy')
    
    # Сохранение
    plt.tight_layout()
    plt.savefig(output_path, 
               dpi=PLOT_SETTINGS['figure']['dpi'],
               bbox_inches='tight',
               pad_inches=0.1,
               facecolor='white')
    plt.close()


def create_plots_2(data, output_path, global_stats):
    """Генерация графиков вдоль линии y = data['Y'][1]"""
    plt.figure(figsize=PLOT_SETTINGS['figure']['size'],
               dpi=PLOT_SETTINGS['figure']['dpi'])

    # Извлечение данных для линии y = data['Y'][1]
    y_index = 1
    x_values = data['X'][y_index, :]
    params_along_y = {
        'vx': data['vx'][y_index, :],
        'vy': data['vy'][y_index, :],
        'p': data['p'][y_index, :],
        'r': data['r'][y_index, :],
        'e': data['e'][y_index, :]
    }

    # Общие настройки
    plt.gca().set_aspect('auto')

    # Создание графиков
    def create_subplot(pos, x_values, y_values, label, title):
        plt.subplot(2, 3, pos)
        plt.plot(x_values, y_values, label=label)
        plt.title(title)
        plt.xlabel('X')
        plt.ylabel(label)
        plt.grid(True)

    # Скорость по X
    create_subplot(1, x_values, params_along_y['vx'], 'Vx', 'X-Velocity')

    # Скорость по Y
    create_subplot(2, x_values, params_along_y['vy'], 'Vy', 'Y-Velocity')

    # Модуль скорости
    plt.subplot(2, 3, 3)
    speed = np.sqrt(params_along_y['vx']**2 + params_along_y['vy']**2)
    plt.plot(x_values, speed, label='Speed')
    plt.title('Speed Magnitude')
    plt.xlabel('X')
    plt.ylabel('Speed')
    plt.grid(True)

    # Давление
    create_subplot(4, x_values, params_along_y['p'], 'Pressure', 'Pressure')

    # Плотность
    create_subplot(5, x_values, params_along_y['r'], 'Density', 'Density')

    # Энергия
    create_subplot(6, x_values, params_along_y['e'], 'Energy', 'Energy')

    # Сохранение
    plt.tight_layout()
    plt.savefig(output_path,
               dpi=PLOT_SETTINGS['figure']['dpi'],
               bbox_inches='tight',
               pad_inches=0.1,
               facecolor='white')
    plt.close()

def create_plots_3(data, output_path, global_stats):
    """Генерация графиков вдоль линии x = data['X'][0][1]"""
    plt.figure(figsize=PLOT_SETTINGS['figure']['size'],
               dpi=PLOT_SETTINGS['figure']['dpi'])

    # Извлечение данных для линии x = data['X'][0][1]
    x_index = 1
    y_values = data['Y'][:, x_index]
    params_along_x = {
        'vx': data['vx'][:, x_index],
        'vy': data['vy'][:, x_index],
        'p': data['p'][:, x_index],
        'r': data['r'][:, x_index],
        'e': data['e'][:, x_index]
    }

    # Общие настройки
    plt.gca().set_aspect('auto')

    # Создание графиков
    def create_subplot(pos, y_values, x_values, label, title, global_min, global_max):
        plt.subplot(2, 3, pos)
        plt.plot(y_values, x_values, label=label)
        plt.ylim(global_min, global_max)  # Фиксация шкалы Y
        plt.title(title)
        plt.xlabel('Y')
        plt.ylabel(label)
        plt.grid(True)

    # Скорость по X
    create_subplot(1, y_values, params_along_x['vx'], 'Vx', 'X-Velocity', global_stats['vx'][0], global_stats['vx'][1])

    # Скорость по Y
    create_subplot(2, y_values, params_along_x['vy'], 'Vy', 'Y-Velocity', global_stats['vy'][0], global_stats['vy'][1])

    # Модуль скорости
    plt.subplot(2, 3, 3)
    speed = np.sqrt(params_along_x['vx']**2 + params_along_x['vy']**2)
    speed_min = np.sqrt(global_stats['vx'][0]**2 + global_stats['vy'][0]**2)
    speed_max = np.sqrt(global_stats['vx'][1]**2 + global_stats['vy'][1]**2)
    plt.plot(y_values, speed, label='Speed')
    plt.ylim(speed_min, speed_max)  # Фиксация шкалы Y
    plt.title('Speed Magnitude')
    plt.xlabel('Y')
    plt.ylabel('Speed')
    plt.grid(True)

    # Давление
    create_subplot(4, y_values, params_along_x['p'], 'Pressure', 'Pressure', global_stats['p'][0], global_stats['p'][1])

    # Плотность
    create_subplot(5, y_values, params_along_x['r'], 'Density', 'Density', global_stats['r'][0], global_stats['r'][1])

    # Энергия
    create_subplot(6, y_values, params_along_x['e'], 'Energy', 'Energy', global_stats['e'][0], global_stats['e'][1])

    # Сохранение
    plt.tight_layout()
    plt.savefig(output_path,
               dpi=PLOT_SETTINGS['figure']['dpi'],
               bbox_inches='tight',
               pad_inches=0.1,
               facecolor='white')
    plt.close()

def create_plots_4(data, output_path, global_stats):
    """Генерация графиков без сетки между клеточками"""
    plt.figure(figsize=PLOT_SETTINGS['figure']['size'], 
              dpi=PLOT_SETTINGS['figure']['dpi'])
    
    # Общие настройки
    plt.gca().set_aspect('equal')
    
    # Создание графиков
    def create_subplot(pos, data, param, label, title):
        plt.subplot(2, 3, pos)
        # Убираем edgecolor и linewidth для исключения сетки
        mesh = plt.pcolormesh(data['X'], data['Y'], data[param],
                            vmin=global_stats[param][0],
                            vmax=global_stats[param][1],
                            cmap='viridis',
                            shading='nearest')
        plt.colorbar(mesh, label=label)
        plt.title(title)
        plt.xlabel('X')
        plt.ylabel('Y')
    # Создание линий уровня
    def create_contour(pos, data, param, label, title, txt_mode = 0):
        fg_color = 'black' if txt_mode == 0 else 'black'
        levels = 40
        plt.subplot(2, 3, pos)
        # Убираем edgecolor и linewidth для исключения сетки
        cont = plt.contour(data['X'], data['Y'], data[param], levels = levels)
        plt.clabel(cont)
        cb = plt.colorbar(cont)
        cb.ax.tick_params(which='both', color=fg_color, labelcolor=fg_color)
        cb.set_label(label, color=fg_color)
        plt.title(title)
        plt.xlabel('X')
        plt.ylabel('Y')
    # Создание линий тока
    def create_streamline(pos, data, param, label, title, txt_mode = 0):
        density = 1.0
        fg_color = 'black' if txt_mode == 0 else 'black'
        plt.subplot(2, 3, pos)
        # Убираем edgecolor и linewidth для исключения сетки
        cont = plt.streamplot(data['X'], data['Y'], data[param[0]], data[param[1]], density = density)
        
        cb = plt.colorbar(cont.lines)
        cb.ax.tick_params(which='both', color=fg_color, labelcolor=fg_color)
        cb.set_label(label, color=fg_color)
        plt.title(title)
        plt.xlabel('X')
        plt.ylabel('Y')
    # Скорость по X
    create_subplot(1, data, 'vx', 'Vx', 'X-Velocity')
    
    # Скорость по Y
    create_subplot(2, data, 'vy', 'Vy', 'Y-Velocity')
    
    # Модуль скорости
    create_streamline(3, data, ['vx', 'vy'], "Velocity", "Velocity streamlines", 1)
    """ plt.subplot(2, 3, 3)
    speed = np.sqrt(data['vx']**2 + data['vy']**2)
    mesh = plt.pcolormesh(data['X'], data['Y'], speed,
                        cmap='viridis',
                        shading='nearest')
    plt.colorbar(mesh, label='Speed')
    plt.title('Speed Magnitude')
    plt.xlabel('X')
    plt.ylabel('Y') """
    
    # Давление
    create_contour(4, data, 'p', 'Pressure', 'Pressure')
    
    # Плотность
    create_contour(5, data, 'r', 'Density', 'Density')
    
    # Энергия
    create_contour(6, data, 'e', 'Energy', 'Energy', 1)
    
    # Сохранение
    plt.tight_layout()
    plt.savefig(output_path.replace('.png', '_no_grid.png'), 
               dpi=PLOT_SETTINGS['figure']['dpi'],
               bbox_inches='tight',
               pad_inches=0.1,
               facecolor='white')
    plt.close()

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # Проверка аргументов
    if len(sys.argv) != 3:
        if rank == 0:
            print("Usage: mpirun -n <procs> python plot_generator.py <input_dir> <output_dir>")
        return
    
    input_dir, output_dir = sys.argv[1], sys.argv[2]
    
    # Создание выходной директории
    if rank == 0:
        os.makedirs(output_dir, exist_ok=True)
    comm.Barrier()
    
    # Получение списка файлов
    if rank == 0:
        files = sorted(glob.glob(os.path.join(input_dir, '*.csv')))
        if not files:
            print(f"No CSV files found in {input_dir}")
            return
    else:
        files = None
    
    files = comm.bcast(files, root=0)
    
    # Распределение файлов
    n_files = len(files)
    chunk = n_files // size
    remainder = n_files % size
    start = rank * chunk + min(rank, remainder)
    end = start + chunk + (1 if rank < remainder else 0)
    local_files = files[start:end]
    
    # Первый проход: сбор статистики
    local_stats = {param: (np.inf, -np.inf) for param in ['p', 'vx', 'vy', 'r', 'e']}
    
    for file_path in local_files:
        try:
            _, stats = read_data(file_path)
            for param in local_stats:
                local_stats[param] = (
                    min(local_stats[param][0], stats[param][0]),
                    max(local_stats[param][1], stats[param][1])
                )
        except Exception as e:
            print(f"[Rank {rank}] Error reading {file_path}: {str(e)}")
    
    # Глобальная статистика
    global_stats = get_global_stats(comm, local_stats)
    
    # Второй проход: генерация графиков
    for file_path in local_files:
        try:
            data, _ = read_data(file_path)
            fname = os.path.basename(file_path).replace('.csv', '.png')
            output_path = os.path.join(output_dir, fname)
            create_plots_4(data, output_path, global_stats)
        except Exception as e:
            print(f"[Rank {rank}] Error processing {file_path}: {str(e)}")
    
    # Завершение
    comm.Barrier()
    if rank == 0:
        print("All plots generated successfully")

if __name__ == "__main__":
    main()