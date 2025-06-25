#!/bin/bash

:> tpr-files.dat
:> pullf-files.dat

ls ../01_firstUSrun/_TPR/umbrella*tpr >> tpr-files.dat
ls ../01_firstUSrun/_XVG/umbrella*pullf.xvg >> pullf-files.dat
ls ../02_secondUSrun/_TPR/umbrella*tpr >> tpr-files.dat
ls ../02_secondUSrun/_XVG/umbrella*pullf.xvg >> pullf-files.dat

gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal -b 0
