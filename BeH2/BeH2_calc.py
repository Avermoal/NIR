from vasp import Vasp
from ase import Atoms, Atom
from ase.io import write
import numpy as np

np.set_printoptions(precision=3, suppress=True)

beh2_bound_length = 1.33#Angstrom exp data

#BeH2 it is linear molecula
atoms = Atoms(
    Atom('H', [0, 0, beh2_bound_length]),
    Atom('Be', [0, 0, 0]),
    Atom('H', [0, 0, -beh2_bound_length]))

L = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
ENCUT = [200, 250, 300, 350, 400, 450, 500]
POTENTIALS = ['LDA', 'PBA', 'PW91', ]

energies = []
ready = True
for xc in POTENTIALS:
    for encut in ENCUT:
        for l in L:
            atoms.set_cell([l, l, l], scale_atoms=False)
            atoms.center()
            calc = Vasp('molecules/beh2-L-{0}'.format(l),
                        encut=encut,
                        xc=xc,
                        atoms=atoms)
            energies.append(atoms.get_potential_energy())

        print(energies)
        calc.stop_if(None in energies)

        import matplotlib.pyplot as plt
        plt.plot(L, energies, 'bo-')
        plt.xlabel('Unit cell length ($\AA$)')
        plt.ylabel('Total energy (eV)')
        plt.savefig('images/beh2-e-v'+xc+encut+'.png')
