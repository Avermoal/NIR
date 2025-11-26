#!/usr/bin/env python

#from jasp import *
from ase.calculators.vasp import Vasp
from ase import Atoms, Atom
from ase.io import write
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

np.set_printoptions(precision=3, suppress=True)

beh2_bound_length = 1.33  # Angstrom exp data

# BeH2 is a linear molecule
atoms = Atoms([
    Atom('H', [0, 0, beh2_bound_length]),
    Atom('Be', [0, 0, 0]),
    Atom('H', [0, 0, -beh2_bound_length])])

L = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
ENCUT = [200, 250, 300, 350, 400, 450, 500]
POTENTIALS = ['LDA', 'PBE', 'PW91']

# Создаем директории для результатов
os.makedirs('molecules', exist_ok=True)
os.makedirs('images', exist_ok=True)

for xc in POTENTIALS:
    for encut in ENCUT:
        energies = []
        for l in L:
            atoms.set_cell([l, l, l], scale_atoms=False)
            atoms.center()
            
            # Создаем уникальное имя директории
            dir_name = 'molecules/beh2-{0}-{1}-{2}'.format(xc, encut, l)
            
            # Очищаем директорию перед расчетом
            if os.path.exists(dir_name):
                shutil.rmtree(dir_name)
            os.makedirs(dir_name)
            
            # Явно создаем входные файлы VASP
            from ase.io.vasp import write_vasp
            write_vasp(os.path.join(dir_name, 'POSCAR'), atoms)
            
            # Создаем калькулятор с явным указанием всех параметров
            calc = Vasp(
                directory=dir_name,
                xc=xc,
                encut=encut,
                prec='Normal',
                ismear=0,
                sigma=0.1,
                ibrion=-1,  # No ionic relaxation
                nsw=0,      # No ionic steps
                lwave=False,
                lcharg=False,
                lreal='Auto',
                nelm=100,
                ediff=1e-6,
                isym=0,
                setups={'Be': '_sv'},  # Используем pseudopotential с semi-core states для Be
                restart=False  # Явно отключаем рестарт
            )
            
            atoms.calc = calc
            
            try:
                energy = atoms.get_potential_energy()
                energies.append(energy)
                print(f"XC: {xc}, ENCUT: {encut}, L: {l}, Energy: {energy}")
            except Exception as e:
                print(f"Error for XC: {xc}, ENCUT: {encut}, L: {l}: {e}")
                energies.append(None)
                # Останавливаем расчеты для этой комбинации, если есть ошибка
                break
        
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
        else:
            print(f"No successful calculations for XC: {xc}, ENCUT: {encut}")
