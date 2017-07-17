# SimCell

BRIEF OVERVIEW
================

This is the README file for SimCell source.


The main repository contains the code of the general geometry engine.

The repository “cells” contains the specific code for segregation and sharpening simulations.




HOW TO USE
================

This code has been designed to work on linux, but you can also run it on MacOS.

Use the makefiles to compile both codes:

cd SimCell

make

cd cells

make



Note: You may need to re-build the links ‘main’ and ‘home’ in cells repository:

cd cells

ln ../ main

ln ./home 




You can then run the code:

./simcell cells-800-ratio_0.5.run




PARAMETERS FILES
================

There are 5 parameters files:
- cell.dat

that contains the description of a cell (number of beads and their geometrical properties). You don’t have to change it.



- cell.param

That contains the parameters of the running simulation (duration, saving step, cell behaviours and motility parameters). You can edit this file to change the simulation duration, turn on/off cell behaviours such as adhesion, repulsion, increase of persistence after contact and modify cell migration parameters (speed, persistence).



- green.cell and red.cell

These two files list the behavioural parameters of each cell populations. You can edit the last two lines to modify the adhesion and repulsion frequencies. 



- cells-800-ratio_0.5.run

This file sets up the cells for the simulation. You can choose their total number (here 00, line GROUP 0 800) and the size of each cell population (number of lines for each model, MODEL 0 and MODEL 1). If you want to change the number of cells, you need to modify the total number and the number of lines for each population (i.e. MODEL).
It allows the code to use the green-cell.model/red-cell.model and cell.dat files by giving it their paths.




CONTACT
================

Please feel free to contact us if you have any question:

William Taylor william.taylor@crick.ac.uk

Anais Khuong anais.khuong@crick.ac.uk

