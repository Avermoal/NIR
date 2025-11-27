#!/usr/bin/env python

#Скрипт вычисляет зависимость потенциальной энергии
#от расстояния между атомами молекулы BeH и записывает
#данные в соответсвующие файлы BeH_PE_(xc)_(encut).txt, где xc --- псевдопотенциал

from ase.calculators.vasp import vasp
from ase.structure import molecule
import numpy as np
import os

atoms = molecule('BeH')

ENCUT = [200, 250, 300, 350, 400, 450, 500]
XC = ['LDA', 'PBE', 'PW91']

#Создаём дерикторию для результатов
os.makedirs('results', exist_ok=True)

for xc in XC:
    for encut in ENCUT:

        energies = []

        for distance in range(0.1, 4, 0.1):

            #Создаём копию atoms для расчётов
            current_atoms = atoms.copy()
            current_atoms.center(vacuum=5)#По-хорошему нужно находить оптимальный размер ячейки в соответсвующих расчётах
            #Устанавливаем расстояние между атомами в молекуле
            current_atoms.set_distance(0, 1, distance)

            try:
                calc = Vasp(
                    'molecules/BeH-L-{0}'.format(distance),
                    encut=encut,
                    xc=xc,
                    atoms=current_atoms,
                    ibrion=-1,#Для отключения имзенения кристаллической структуры т.е. без релаксации
                    nsw=0#Нуль ionic steps
                )

                current_atoms.calc = calc

                energy = current_atoms.get_potential_energy()
                energies.append(energy)

            except Exception as e:
                print(f"Error for XC: {xc}, ENCUT: {encut}, L: {l}: {e}")
                energies.append(None)

        #Пишем только успешные расчёты
        if any(e is not None for e in energies):
            # Фильтруем None значения для построения графика
            valid_L = [L[i] for i in range(len(energies)) if energies[i] is not None]
            valid_energies = [e for e in energies if e is not None]
            
            # Сохраняем данные в текстовый файл
            with open(f'results/beh-e-v-{xc}-{encut}.txt', 'w') as f:
                for i, l in enumerate(valid_L):
                    f.write(f"{l} {valid_energies[i]:.6f}\n")
        else:
            print(f"No successful calculations for XC: {xc}, ENCUT: {encut}")

print("Расчет завершен!")

