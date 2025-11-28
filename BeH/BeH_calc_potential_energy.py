#!/usr/bin/env python

#Скрипт вычисляет зависимость потенциальной энергии
#от расстояния между атомами молекулы BeH и записывает
#данные в соответсвующие файлы BeH_PE_(xc)_(encut).txt, где xc --- псевдопотенциал

from ase.calculators.vasp import Vasp
from ase.structure import molecule
import numpy as np
import os

atoms = molecule('BeH')

ENCUT = [200, 250, 300, 350, 400, 450, 500]
XC = ['LDA', 'PBE']

#Создаём дерикторию для результатов
os.makedirs('results', exist_ok=True)

for xc in XC:
    for encut in ENCUT:

        energies = []
        dist = []

        for distance in np.arange(0.1, 4, 0.1):
            dist.append(distance)
            #Создаём копию atoms для расчётов
            current_atoms = atoms.copy()
            current_atoms.center(vacuum=5)#По-хорошему нужно находить оптимальный размер ячейки в соответсвующих расчётах
            #Устанавливаем расстояние между атомами в молекуле
            current_atoms.set_distance(0, 1, distance)

                # Создаем уникальную директорию для каждого расчета
            calc_dir = f'molecules/BeH-{xc}-{encut}-{distance:.2f}'
            
            # Удаляем директорию, если она существует, чтобы начать "чистый" расчет
            if os.path.exists(calc_dir):
                shutil.rmtree(calc_dir)
            
            try:

                               # Создаем калькулятор с улучшенными настройками для молекул
                calc = Vasp(
                    encut=encut,
                    xc=xc,
                    ibrion=-1,     # Без релаксации ионов
                    nsw=0,          # Без шагов ионной релаксации
                    istart=0,       # Новый расчет
                    icharg=2,       # Автоматическое определение начальной зарядовой плотности
                    lwave=False,    # Не записывать волновые функции
                    lcharg=False,   # Не записывать зарядовую плотность
                    prec='Normal',  # Точность расчета
                    ismear=0,       # Гауссово сглаживание (подходит для молекул)
                    sigma=0.01,     # Малая ширина сглаживания для молекул
                    ediff=1e-6,     # Точность по энергии
                    nelm=100,       # Максимальное число электронных шагов
                    nelmin=4,       # Минимальное число электронных шагов
                    algo='Normal',  # Алгоритм электронной минимизации
                    lreal='Auto',   # Автоматический выбор реального/обратного пространства
                    isym=0,         # Без симметрии (для молекул)
                    lorbit=11,      # Вывод информации о проекциях
                    npar=1,         # Последовательный расчет (без параллелизации)
                    kpts=(1, 1, 1)  # Только гамма-точка для молекул
                )
                
                current_atoms.calc = calc
                # Устанавливаем директорию отдельно
                calc.directory = calc_dir

                current_atoms.calc = calc

                energy = current_atoms.get_potential_energy()
                energies.append(energy)

            except Exception as e:
                print(f"Error for XC: {xc}, ENCUT: {encut}, DISTANCE: {distance}: {e}")
                energies.append(None)

        #Пишем только успешные расчёты
        if any(e is not None for e in energies):
            # Фильтруем None значения для построения графика
            valid_dist = [dist[i] for i in range(len(energies)) if energies[i] is not None]
            valid_energies = [e for e in energies if e is not None]
            
            # Сохраняем данные в текстовый файл
            with open(f'results/beh-e-v-{xc}-{encut}.txt', 'w') as f:
                for i, l in enumerate(valid_dist):
                    f.write(f"{l} {valid_energies[i]:.2f}\n")
        else:
            print(f"No successful calculations for XC: {xc}, ENCUT: {encut}")

print("Расчет завершен!")

