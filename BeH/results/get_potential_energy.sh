#!bin/bash

#Создаёт графики зависимости потенциальной энергии от расстояния между атомами в молекуле BeH ENCUT для псевдопотенциалов LDA и PBE

echo "set terminal pdf
set output 'BeH_potenergy_LDA.pdf'
set xlabel 'Interatomic distance, Angstrom'
set ylabel 'Potential energy, eV'
set key right top
set key font 'Arial, 7'
plot 'beh-e-v-LDA-200.txt' w lp title 'ENCUT = 200', 'beh-e-v-LDA-250.txt' w lp title 'ENCUT = 250', 'beh-e-v-LDA-300.txt' w lp title 'ENCUT = 300', 'beh-e-v-LDA-350.txt' w lp title 'ENCUT = 350', 'beh-e-v-LDA-400.txt' w lp title 'ENCUT = 400', 'beh-e-v-LDA-450.txt' w lp title 'ENCUT = 450', 'beh-e-v-LDA-500.txt' w lp title 'ENCUT = 500'," | gnuplot

echo "set terminal pdf
set output 'BeH_potenergy_PBE.pdf'
set xlabel 'Interatomic distance, Angstrom'
set ylabel 'Potential energy, eV'
set key right top
set key font 'Arial, 7'
plot 'beh-e-v-PBE-200.txt' w lp title 'ENCUT = 200', 'beh-e-v-PBE-250.txt' w lp title 'ENCUT = 250', 'beh-e-v-PBE-300.txt' w lp title 'ENCUT = 300', 'beh-e-v-PBE-350.txt' w lp title 'ENCUT = 350', 'beh-e-v-PBE-400.txt' w lp title 'ENCUT = 400', 'beh-e-v-PBE-450.txt' w lp title 'ENCUT = 450', 'beh-e-v-PBE-500.txt' w lp title 'ENCUT = 500'," | gnuplot
