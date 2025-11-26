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
    Atom('H', [0, 0, -beh2_bound_length])])

# Уменьшим размеры сетки для тестирования
L = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
ENCUT = [200, 250, 300, 350, 400, 450, 500]
POTENTIALS = ['LDA', 'PBE', 'PW91']

# Создаем директории для результатов
os.makedirs('images', exist_ok=True)

for xc in POTENTIALS:
    for encut in ENCUT:
        energies = []
        for l in L:
            atoms.set_cell([l, l, l], scale_atoms=False)
            atoms.center()
            
            # Создаем временную директорию
            with tempfile.TemporaryDirectory() as tmp_dir:
                try:
                    # Минимальная конфигурация - только самые необходимые параметры
                    calc = Vasp(
                        tmp_dir,
                        xc=xc,
                        encut=encut,
                        restart = None
                    )
                    
                    atoms.calc = calc
                    
                    energy = atoms.get_potential_energy()
                    energies.append(energy)
                    print(f"XC: {xc}, ENCUT: {encut}, L: {l}, Energy: {energy}")
                except Exception as e:
                    print(f"Error for XC: {xc}, ENCUT: {encut}, L: {l}: {e}")
                    energies.append(None)
        
        # Строим график только если есть успешные расчеты
        if any(e is not None for e in energies):
            # Фильтруем None значения для построения графика
            valid_L = [L[i] for i in range(len(energies)) if energies[i] is not None]
            valid_energies = [e for e in energies if e is not None]
            
            plt.figure()
            plt.plot(valid_L, valid_energies, 'bo-')
            plt.xlabel('Unit cell length ($\AA$)')
            plt.ylabel('Total energy (eV)')
            plt.title(f'BeH2 - XC: {xc}, ENCUT: {encut}')
            plt.savefig(f'images/beh2-e-v-{xc}-{encut}.png')
            plt.close()
            
            # Сохраняем данные в текстовый файл
            with open(f'images/beh2-e-v-{xc}-{encut}.txt', 'w') as f:
                for i, l in enumerate(valid_L):
                    f.write(f"{l} {valid_energies[i]}\n")
        else:
            print(f"No successful calculations for XC: {xc}, ENCUT: {encut}")
