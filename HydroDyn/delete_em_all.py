import os
import shutil

def clear_folder_contents(folder_path):
    """
    Удаляет всё содержимое указанной папки, не удаляя саму папку.
    """
    if os.path.exists(folder_path):
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            try:
                if os.path.isfile(item_path) or os.path.islink(item_path):
                    os.unlink(item_path)  # Удалить файл или ссылку
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)  # Удалить папку
            except Exception as e:
                print(f"Не удалось удалить {item_path}: {e}")
        print(f"Содержимое папки {folder_path} успешно удалено.")
    else:
        print(f"Папка {folder_path} не существует.")

# Запрос номеров папок
input_numbers = input("Введите номера папок через запятую (например, 1,2,3): ").strip()

# Преобразование ввода в список номеров
try:
    folder_numbers = [int(num.strip()) for num in input_numbers.split(',') if num.strip().isdigit()]
except ValueError:
    print("Неверный формат ввода. Укажите номера папок через запятую.")
    exit(1)

# Обработка каждой папки
for num in folder_numbers:
    folder_name = f"results{num}"
    clear_folder_contents(folder_name)
