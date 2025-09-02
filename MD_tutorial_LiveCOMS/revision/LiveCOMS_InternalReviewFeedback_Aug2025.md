## Tristan Feedback
 - [X] P1: „The Collection of all Forces acting on all atoms is called the Forcefield and forms the primary assumption of each MD simulation.” ForceField is the Assumption that PotEner can be approximated by analytic functions with empirical parameters. Therefore the quoted text is misleading. 
   > [!NOTE] No changes were applied
   > The collection of all terms and params is the collection of the potentials 
   > and the collection of the potentials is the collection of the forcces (F=dU/dr)
   > The Introduction is supposed to be brief and is not suitable for explaining stuff in detail
 - [X] P2: Correction: Spelling error: effects -> affect;
   > [!DONE] Correction Conducted
 - [X] P2: “especially” seems to be a bit too much. “These limitations affect the physical realism of MD simulations, with the accurate treatment of protonation state fluctuations being a prominent problem in the study of biological systems. …”
   > [!NOTE] No changes were applied
   > Protonation is the most frequent reaction that takes place
   > important for all processes, also thos requiring chemistry
   > e.g. catalytic reactions often feature prot/deprot of solvent molecules before catalysis
 - [X] P3: GROMACS 2026?
   > [!NOTE] No changes were applied
   > yes github dev version is 2026
 - [X] P7: no numbering on Fig 5
   > [!DONE] Revision Incorporated
   > Merged the subfigures into one figure and changed the caption
 - [X] P7: First correction: Otherwise, A trajectory … -> Otherwise, a trajectory 
   > [!DONE] Correction Conducted
 - [X] P8: “Assuming, N particles the motions governing particle movement [...].” Difficult to understand the meaning. Maybe: “motions” -> laws or equations
   > [!DONE] Revision Incorporated
   > The sentence was indeed convoluted
   > changed to: "Assuming $N$ particles, the particle movement depends on the forces $F$ on each particle and the potential energy function [...]"
 - [X] P8: Next sentence: “The potential energy function is also called the forcefield …” - This might lead to confusion.
   > [!DONE] Revision Incorporated
   > The two are equivalent terms as described above. I changed it to 'total potential energy function'
 - [X] P8 (Also P16,23,27): GMX help link: if “documentation/2018” is switched to „current“ it should Point to the newest gmx documentation.
   > [!DONE] Revision incorporated
   > cmdline changed to https://manual.gromacs.org/current/user-guide/cmdline.html
   > manual changed to https://manual.gromacs.org/current/reference-manual/index.html
   > mdpoptions changed to https://manual.gromacs.org/current/user-guide/mdp-options.html
   > gmx mdrum changed to https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html
 - [X] P9: Table, First line: Misses a Point After the sentence defining .tpr
   > [!DONE] Revision incorporated
 - [X] P33: Left side:, under figure: separation line seems to be in code format 
   > [!WARNING] Unclean Fix
   > No idea what caused it but placing [H] behind the \begin figure statement fixed it (no line at all)
   > might be a good idea to try their new cls / bst file to see if that fixes it
   > I have a feeling the bug might stem from their setup
 - [X] P25: MBAR method: Do we actually use MBAR here, as gmx bar i think only uses BAR.
   > [!NOTE] No changes were applied
   > BAR is a special case of MBAR so technically by stating MBAR, BAR is included

---

## Martin Feedback 
 - [X] AMBER, NAMD etc. are getting mentioned, so they need to be cited.
   > [!DONE] Revision Incorporated
   > citations were at the end of the sentence
   > now they are moved to the specific simulation software for clearer citation
 - [X] Modeling „histidine“ can be achieved with MD packages that implemented constant pH simulations. Not strictly necessary to mention.
   > [!NOTE] No changes were applied
   > I feel like this is explaining too much, a feature/phenomenon that is not important for the tutorial itself
 - [ ] Use martins response letter
---

## Shuyus Feedback
- [X] Figure 44 is a little too small in the pdf. make it single column
  > [!DONE] Revision Incorporated
  > Image is now a single column picture
- [X] The link medium link to the smiles tutorial does not work. Instead I can access the blog via [...]
  > [!DONE] Revision Incorporated
  > new link https://luis-vollmers.medium.com/tutorial-to-smiles-and-canonical-smiles-explained-with-examples-fbc8a46ca29f
- [X] right figure of Figure 48 is the zoomed-in view of the left instead of another docking result. include this in the caption as well.
  > [!DONE] Revision Incorporated
- [X] On page 54: "Therefore, it is recommended not to use mv rather than rm"  --> I assume the "not" here is a typo?
  > [!DONE] Revision Incorporated
- [X] Figure 54-57 mixing color schemes. distinguish different repeats but always use light blue for compound 0156 and orange for 1258.
  > [!DONE] Revision Incorporated
  > I used the green/purple combo and used different line/marker styles to distinguish runs. Looks neat now imo.

---

