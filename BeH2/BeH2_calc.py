#!/usr/bin/env python

from ase.calculators.vasp import Vasp
from ase import Atoms, Atom
import numpy as np
import matplotlib.pyplot as plt
import os
import tempfile

np.set_printoptions(precision=3, suppress=True)

beh2_bound_length = 1.33  # Angstrom exp data

# BeH2 is a linear molecule
atoms = Atoms([
    Atom('H', [0, 0, beh2_bound_length]),
    Atom('Be', [0, 0, 0]),
    Atom('H', [0, 0, -beh2_bound_length])],
    pbc=True)

L = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
ENCUT = [200, 250, 300, 350, 400, 450, 500]
POTENTIALS = ['LDA', 'PBE']

# Создаем директории для результатов
os.makedirs('images', exist_ok=True)

print('atom symbol')
print('===========')
for i, atom in enumerate(atoms):
    print('{0:2d} {1:3s}'.format(i, atom.symbol))

for xc in POTENTIALS:
    for encut in ENCUT:
        energies = []
        angles = []
        for l in L:
            # Создаем копию атомов для каждого расчета
            current_atoms = atoms.copy()
            current_atoms.set_cell([l, l, l], scale_atoms=False)
            current_atoms.center()
            
            # Создаем временную директорию
            with tempfile.TemporaryDirectory() as tmp_dir:
                try:
                    # Базовая конфигурация VASP
                    calc = Vasp(
                        xc=xc,
                        encut=encut,
                        istart=0,  # начинать с чистого листа
                        icharg=2,  # атомные заряды
                        prec='Normal',
                        ismear=0,  # Гауссово размазывание
                        sigma=0.1,  # ширина размазывания
                        ediff=1e-4,  # точность по энергии
                        nsw=0,  # без релаксации ионов
                        ibrion=-1,  # без релаксации
                        isif=2,  # рассчитать напряжения
                        lwave=False,  # не писать файлы волновых функций
                        lcharg=False,  # не писать файлы зарядов
                        lreal='Auto'
                    )

                    calc.directory=tmp_dir
                    
                    current_atoms.calc = calc
                    
                    energy = current_atoms.get_total_energy()
                    energies.append(energy)
                    print(f"XC: {xc}, ENCUT: {encut}, L: {l}, Energy: {energy:.6f}")

                    # Угл связи
                    da = atoms.get_dihedral([1, 0, 2]) * 180. / np.pi
                    angles.append(da)
                    
                except Exception as e:
                    print(f"Error for XC: {xc}, ENCUT: {encut}, L: {l}: {e}")
                    energies.append(None)
        
        # Строим график только если есть успешные расчеты
        if any(e is not None for e in energies):
            # Фильтруем None значения для построения графика
            valid_L = [L[i] for i in range(len(energies)) if energies[i] is not None]
            valid_energies = [e for e in energies if e is not None]

            # Сохраняем данные в текстовый файл
            with open(f'images/beh2-e-v-{xc}-{encut}.txt', 'w') as f:
                for i, l in enumerate(valid_L):
                    f.write(f"{l} {valid_energies[i]:.6f} {angles[i]:.6f}\n")
        else:
            print(f"No successful calculations for XC: {xc}, ENCUT: {encut}")

print("Расчет завершен!")
