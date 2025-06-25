#!/bin/bash

system_preparation()
{
  # initialise files data files for wham
  :> tpr-files.dat
  :> pullf-files.dat
  # create topologies and coordinates for GROMACS from PDB file
  printf "13\n1\n2\n0\n2\n0\n2\n0\n2\n0\n2\n0\n" | gmx pdb2gmx -f 2BEG_model1_capped.pdb -ignh -ter -o complex.gro
  # this line is necessary for the pulling simulation
  printf "#ifdef POSRES_B\n#include \"posre_Protein_chain_B.itp\"\n#endif\n" >> topol_Protein_chain_B.itp
  # place the protofibril at the correct location and create the box, values from online tutorial.
  gmx editconf -f complex.gro -o newbox.gro -center 3.280 2.181 2.4775 -box 6.560 4.362 12
  # solvate the newly created box
  gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
  # ion generation in GROMACS requires a TPR file, therefore use grompp then genion...
  gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
  printf "13\n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1
}

box_equilibration()
{
  # run minimisation and equilibration
  gmx grompp \
      -f em.mdp \
      -c solv_ions.gro \
      -p topol.top \
      -o em.tpr \
      -maxwarn 1
  gmx mdrun -v -deffnm em
  gmx grompp \
      -f npt.mdp \
      -c em.gro \
      -p topol.top \
-r em.gro \
      -o npt.tpr \
      -maxwarn 3
  gmx mdrun -deffnm npt
  # make an index that contains chain a and chain b explicitly
  gmx make_ndx -f npt.gro<<EOF
r 1-27
r 28-54
q
EOF
  sed -i 's/r_1-27/Chain_A/g' index.ndx
  sed -i 's/r_28-54/Chain_B/g' index.ndx
}

pull_simulation()
{
  # Run the pulling simulation
  gmx grompp \
      -f md_pull.mdp \
      -c npt.gro \
      -p topol.top \
      -r npt.gro \
      -n index.ndx \
      -t npt.cpt \
      -o pull.tpr \
      -maxwarn 3
  gmx mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg
  gmx pairdist -f pull.xtc -s pull.tpr -n index.ndx -sel "com of group Chain_A" -ref "com of group Chain_B" -o dist.xvg
  sed '/^\@/d' dist.xvg |\
  sed '/^\#/d' |\
  awk '$2 > i+0.2{a[k++]=$1;i=i+0.2}END{for(j=0;j<501;j++)if(a[j] != 0.000 && j > 0 || j < 1 )printf("%7.3f%1s",a[j],"\n")}' > important_frames.txt       # write out every frame that has COM dist 0.2 higher than the prior one
  # if another spacing is wanted the i+0.2 can be changed into i+spacing
}

umbrella_sampling_simulations()
{
  # generate the coordinate files from the important frames
  for frame in $(cat important_frames.txt)
  do
          printf "0\n" | gmx trjconv -s pull.tpr -f pull.xtc -o conf${frame%%.*}.gro -b $frame -e $frame
  done
  for frame in $(cat important_frames.txt)
  do
          i=${frame%%.*}
          echo "running for frame $i..."
          gmx grompp \
              -f npt_umbrella.mdp \
              -c conf$i.gro \
              -p topol.top \
              -r conf$i.gro \
              -n index.ndx \
              -o npt$i.tpr \
              -maxwarn 3
          gmx mdrun -v -deffnm npt$i
          gmx grompp \
              -f md_umbrella.mdp \
              -c npt$i.gro \
              -t npt$i.cpt \
              -p topol.top \
              -r npt$i.gro \
              -n index.ndx \
              -o umbrella$i.tpr \
              -maxwarn 1
          gmx mdrun -v -deffnm umbrella$i
          echo "umbrella$i.tpr" >> tpr-files.dat
          echo "umbrella${i}_pullf.xvg" >> pullf-files.dat
          echo "done"
  done

  gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal -b 0
}

additional_windows_simulations()
{
extraframe=( 50 147 186 195 257 302 340 393 ) # already added: ...
if [[ -z "$extraframe" ]]; then
   echo "no frames specified to run an umbrella simulation for." ; exit 1
else
        for i in "${extraframe[@]}"
        do
                echo "running for frame $i"
                printf "0\n" | gmx trjconv -s pull.tpr -f pull.xtc -o conf$i.gro -b $i -e $i
                gmx grompp \
                    -f npt_umbrella.mdp \
                    -c conf$i.gro \
                    -p topol.top \
                    -r conf$i.gro \
                    -n index.ndx \
                    -o npt$i.tpr \
                    -maxwarn 3
                gmx mdrun -v -deffnm npt$i
                gmx grompp \
                    -f md_umbrella.mdp \
                    -c npt$i.gro \
                    -t npt$i.cpt \
                    -p topol.top \
                    -r npt$i.gro \
                    -n index.ndx \
                    -o umbrella$i.tpr \
                    -maxwarn 1
                gmx mdrun -v -deffnm umbrella$i
                echo "umbrella$i.tpr" >> tpr-files.dat
                echo "umbrella${i}_pullf.xvg" >> pullf-files.dat
                echo "done"
        done
fi
gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal -b 0
}

cleanup()
{
  mkdir -p _CPT  _EDR  _GRO  _LOG  _TPR  _trash  _TRR  _XTC  _XVG
  mv -v *.cpt _CPT
  mv -v *.edr _EDR
  mv -v *.gro _GRO
  mv -v *.log _LOG
  mv -v *.tpr _TPR
  mv -v *.trr _TRR
  mv -v *.xtc _XTC
  mv -v *.xvg _XVG
  mv -v \#* _trash
}

timings()
{
  echo "total/s per.CPU/s Acc" | awk '{printf "%12s %11s %6s\n", $1,$2,$3}' > timings.txt
  for log in $(ls _LOG/*)
  do
    grep 'Time:' $log | awk '{printf "%12.3f %11.3f %6.1f\n", $2,$3,$4}' >> timings.txt
  done
  sum_total=$(awk 'NR > 1{print $1}' timings.txt | paste -sd+ | bc)
  sum_pCPU=$(awk 'NR > 1{print $2}' timings.txt | paste -sd+ | bc)
  echo "$sum_total $sum_pCPU" | awk '{printf "%12.3f %11.3f",$1,$2}' >> timings.txt
}

echo "Make sure to uncomment functions before executing."
#system_preparation
#box_equilibration
#pull_simulation
#umbrella_sampling_simulations
#additional_windows_simulations
#cleanup
#timings
