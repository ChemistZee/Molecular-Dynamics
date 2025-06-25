#!/bin/bash
# Gromacs Formula for calculating Isothermal Compressibility k_T
#
#        RMSD(Vol) * RMSD(Vol)                                                          
# k_T = -----------------------                                                                     
#           k_B * T * <Vol>                                                               

solve_first_task()
{
  # make summary file
  gmx --version > summary.txt
  printf "\n\nExec 01 - Barostat Comparison" >> summary.txt
  # create system: solvate -> insert -> pdb2gmx
  gmx solvate -cs spc216 -box 2.504 2.504 2.504 -maxsol 1 -o tip3p.gro
  gmx insert-molecules -ci tip3p.gro -box 2.504 2.504 2.504 -nmol 523 -try 10000 -o choice_tauP.gro
  printf "6\n" | gmx pdb2gmx -f choice_tauP.gro -p choice_tauP.top -water tip3p -o _trash_pdb2gmx.gro 
  # EM
  gmx grompp \
      -f em.mdp \
      -c choice_tauP.gro \
      -p choice_tauP.top \
      -o choice_tauP_em.tpr \
      -maxwarn 1
  gmx mdrun -s choice_tauP_em.tpr -deffnm choice_tauP_em
  # NVT
  gmx grompp \
      -f nvt.mdp \
      -c choice_tauP_em.gro \
      -p choice_tauP.top \
      -o choice_tauP_nvt.tpr \
      -maxwarn 1
  gmx mdrun -s choice_tauP_nvt.tpr -deffnm choice_tauP_nvt
  for bstat in Berendsen Parrinello-Rahman
  do
    printf "\n$bstat\n" >> summary.txt
    sed "s/BAROSTAT/$bstat/g" choice_npt.mdp > choice_npt_${bstat:0:3}.mdp
    for taup in 0.05 0.1 0.5 1 5 10
    do
      printf "tauP:    $taup\n" >> summary.txt
      printf "RMSD(Vol):   " >> summary.txt
  
      sed "s/TAUP/$taup/g" choice_npt_${bstat:0:3}.mdp > choice_npt_${bstat:0:3}_${taup}.mdp
      if (( $(echo "$taup > 0.1" | bc -l) )); then sed -i '/^nstpcouple.*/d' choice_npt_Par_${taup}.mdp; fi
      sed -i '/^nstpcouple.*/d' choice_npt_Ber_*.mdp
      gmx grompp \
          -f choice_npt_${bstat:0:3}_${taup}.mdp \
          -c choice_tauP_nvt.gro -p choice_tauP.top \
          -o choice_npt_${bstat:0:3}_${taup}.tpr \
          -maxwarn 2
      gmx mdrun -v -s choice_npt_${bstat:0:3}_${taup}.tpr -deffnm choice_npt_${bstat:0:3}_${taup}
      printf "Volume\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o vol_choice_npt_${bstat:0:3}_${taup}.xvg | grep 'nm^3' | awk '{print $4}' >> summary.txt
      printf "Pressure\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o pres_choice_npt_${bstat:0:3}_${taup}.xvg
      printf "Temperature\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o temp_choice_npt_${bstat:0:3}_${taup}.xvg
      printf "Total-Energy\n" |  gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o etot_choice_npt_${bstat:0:3}_${taup}.xvg
      printf "Volume\nTemperature" | gmx energy -f choice_npt_${bstat:0:3}_${taup}.edr -o _trash_temp.xvg -fluct_props -nmol 523 | egrep -o 'Kappa.*' >> summary.txt
    done
  done
}

solve_second_task()
{
  printf "\n\nExec 02 - Watermodel Comparison\n" >> summary.txt
  # create system: solvate -> insert -> pdb2gmx
  gmx solvate -cs tip4p -box 2.504 2.504 2.504 -maxsol 1 -o tip4p.gro
  gmx solvate -cs spc216 -box 2.504 2.504 2.504 -maxsol 1 -o spce.gro
  gmx insert-molecules -ci tip4p.gro -box 2.504 2.504 2.504 -nmol 523 -try 10000 -o tip4pbox.gro
  gmx insert-molecules -ci spce.gro -box 2.504 2.504 2.504 -nmol 523 -try 10000 -o spcebox.gro
  printf "6\n" | gmx pdb2gmx -f tip4pbox.gro -p tip4pbox.top -water tip4p -o _trash_pdb2gmx.gro 
  printf "6\n" | gmx pdb2gmx -f spcebox.gro -p spcebox.top -water spce -o _trash_pdb2gmx.gro 
  # EM
  gmx grompp -f em.mdp -c spcebox.gro -p spcebox.top -o spcebox_em.tpr
  gmx grompp -f em.mdp -c tip4pbox.gro -p tip4pbox.top -o tip4pbox_em.tpr
  gmx mdrun -s spcebox_em.tpr -deffnm spcebox_em
  gmx mdrun -s tip4pbox_em.tpr -deffnm tip4pbox_em
  # NVT
  gmx grompp -f nvt.mdp -c spcebox_em.gro -p spcebox.top -o spcebox_nvt.tpr -maxwarn 1
  gmx grompp -f nvt.mdp -c tip4pbox_em.gro -p tip4pbox.top -o tip4pbox_nvt.tpr -maxwarn 1
  gmx mdrun -s spcebox_nvt.tpr -deffnm spcebox_nvt
  gmx mdrun -s tip4pbox_nvt.tpr -deffnm tip4pbox_nvt
  # NPT
  gmx grompp -f waters_npt.mdp -c spcebox_nvt.gro -p spcebox.top -o spcebox_npt.tpr -maxwarn 2
  gmx grompp -f waters_npt.mdp -c tip4pbox_nvt.gro -p tip4pbox.top -o tip4pbox_npt.tpr -maxwarn 2
  gmx mdrun -s spcebox_npt.tpr -deffnm spcebox_npt
  gmx mdrun -s tip4pbox_npt.tpr -deffnm tip4pbox_npt
  # PRE
  gmx grompp -f waters_pre.mdp -c spcebox_npt.gro -p spcebox.top -o spcebox_pre.tpr
  gmx grompp -f waters_pre.mdp -c tip4pbox_npt.gro -p tip4pbox.top -o tip4pbox_pre.tpr
  gmx mdrun -s spcebox_pre.tpr -deffnm spcebox_pre
  gmx mdrun -s tip4pbox_pre.tpr -deffnm tip4pbox_pre
  # PROD
  gmx grompp -f waters_prod.mdp -c spcebox_pre.gro -p spcebox.top -o spcebox_prod.tpr
  gmx grompp -f waters_prod.mdp -c tip4pbox_pre.gro -p tip4pbox.top -o tip4pbox_prod.tpr
  gmx mdrun -v -s spcebox_prod.tpr -deffnm spcebox_prod
  gmx mdrun -v -s tip4pbox_prod.tpr -deffnm tip4pbox_prod
  # DATA ANALYSIS
  printf "TIP4P " >> summary.txt
  printf "Volume\nTemperature\n" | gmx energy -f tip4pbox_prod.edr -o tip4pbox_prod.xvg -fluct_props -nmol 523 | egrep -o 'Kappa.*' >> summary.txt
  printf "SPC/E " >> summary.txt
  printf "Volume\nTemperature\n" | gmx energy -f spcebox_prod.edr -o spcebox_prod.xvg -fluct_props -nmol 523 | egrep -o 'Kappa.*' >> summary.txt
}

solve_third_task()
{
  # SOLUTION: The ensemble that is simulated is no longer the NPT ensemble but the NPH ensemble. 
  # SOLUTION: Therefore the energy values drift away from the ensemble averages of the NPT ensemble.
  # SOLUTION: A way to test this hypothesis is to calculate properties related to the NPH ensemble like the Joule-Thompson Coefficient.
  printf "\n\nExec 03 - Breaking the System\n" >> summary.txt
  # NPH
  gmx grompp -f standard.mdp -c spcebox_nvt.gro -p spcebox.top -o standard.tpr -maxwarn 2
  gmx mdrun -v -s standard.tpr -deffnm standard
  gmx grompp -f non_standard.mdp -c spcebox_nvt.gro -p spcebox.top -o non_standard.tpr -maxwarn 3
  gmx mdrun -v -s non_standard.tpr -deffnm non_standard
  printf "\nstandard\n" >> summary.txt
  printf "Density\n" | gmx energy -f standard.edr -o dens_standard.xvg | egrep '(kg/m\^3)' >> summary.txt
  printf "Enthalpy\n" | gmx energy -f standard.edr -o enth_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "Potential\n" | gmx energy -f standard.edr -o epot_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "Kinetic-En.\n" | gmx energy -f standard.edr -o ekin_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "Total-Energy\n" | gmx energy -f standard.edr -o etot_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "\nnon_standard\n" >> summary.txt
  printf "Density\n" | gmx energy -f non_standard.edr -o dens_non_standard.xvg | egrep '(kg/m\^3)' >> summary.txt
  printf "Enthalpy\n" | gmx energy -f non_standard.edr -o enth_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "Potential\n" | gmx energy -f non_standard.edr -o epot_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "Kinetic-En.\n" | gmx energy -f non_standard.edr -o ekin_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
  printf "Total-Energy\n" | gmx energy -f non_standard.edr -o etot_non_standard.xvg | egrep '(kJ/mol)' >> summary.txt
}

cleanup()
{
  mkdir -p _trash _output
  mv -v temp.xvg *.trr *.cpt *.edr *.log *.tpr *.xtc mdout.mdp _output
  mv -v step*.pdb \#* mdout.mdp _trash
}


echo "Make sure to uncomment functions before executing."
#solve_first_task
#solve_second_task
#solve_third_task
cleanup
