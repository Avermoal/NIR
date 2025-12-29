#!bin/bash

#Создаёт графики зависимости полной энергии от размеров ячейки BeH2 ENCUT для псевдопотенциалов LDA и PBE

echo "set terminal pdf
set output 'BeH2_totenergy_LDA.pdf'
set xlabel 'Unit cell, Angstrom'
set ylabel 'Total energy, eV'
set key center top
set key font 'Arial, 7'
plot 'beh2-e-v-LDA-200.txt' w lp title 'ENCUT = 200', 'beh2-e-v-LDA-250.txt' w lp title 'ENCUT = 250', 'beh2-e-v-LDA-300.txt' w lp title 'ENCUT = 300', 'beh2-e-v-LDA-350.txt' w lp title 'ENCUT = 350', 'beh2-e-v-LDA-400.txt' w lp title 'ENCUT = 400', 'beh2-e-v-LDA-450.txt' w lp title 'ENCUT = 450', 'beh2-e-v-LDA-500.txt' w lp title 'ENCUT = 500'," | gnuplot

echo "set terminal pdf
set output 'BeH2_totenergy_PBE.pdf'
set xlabel 'Unit cell, Angstrom'
set ylabel 'Total energy, eV'
set key right top
set key font 'Arial, 7'
plot 'beh2-e-v-PBE-200.txt' w lp title 'ENCUT = 200', 'beh2-e-v-PBE-250.txt' w lp title 'ENCUT = 250', 'beh2-e-v-PBE-300.txt' w lp title 'ENCUT = 300', 'beh2-e-v-PBE-350.txt' w lp title 'ENCUT = 350', 'beh2-e-v-PBE-400.txt' w lp title 'ENCUT = 400', 'beh2-e-v-PBE-450.txt' w lp title 'ENCUT = 450', 'beh2-e-v-PBE-500.txt' w lp title 'ENCUT = 500'," | gnuplot
