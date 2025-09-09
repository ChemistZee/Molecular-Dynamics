# MD Tutorial: From Basics to Application

This directory contains the tutorial/training hybrid article "Molecular Dynamics: From Basics to Application" (for authors, see the list below). This document is designed to teach experienced undergraduate or graduate students the basics of molecular dynamics simulations and enable them to tackle their computational biophysics projects. Furthermore, this repository allows any user to contact the authors, report bugs, and participate in improving this article.

## The Vision

Molecular dynamics simulations are a powerful computational tool for understanding the physics underlying biochemical processes. The central vision of this tutorial is to empower anyone with the appropriate scientific background (or sheer will power) to learn these skills independently. At the same time, the tutorial is designed to be valuable for experienced researchers. Many of the authors are themselves seasoned scientists, yet still gained new insights while compiling this work.
A major part of the vision is to make the tutorial didactically meaningful through deliberate design choices:
 - Information is introduced when it becomes relevant, with milestone research publications cited for further reading.
 - Learning objectives are not handed over directly but are meant to be discovered proactively through exercises.
 - Simulation algorithms are not presented as an end in itself but are explored as part of physical problem solving.
 - The tutorial teaches not just which algorithmic options to choose, but how to verify those choices self-reliantly in the context of molecular dynamics.
 - Every exercise subdirectoriy contains the exercise material (```student/```) and the solution material (```teacher/```).
 - ```student/``` folders only contain the starter files required to complete the exercise tasks self-reliantly.
 - ```teacher/``` folders contain all solutions, intermediate files, and executables that automatically generate the solutions from the starter files.
Finally, the vision is to encourage university lecturers to adopt this tutorial in full or in part for their classes. By fostering positive feedback loops among lecturers, students, and independent users, this project can grow as a living publication that continually improves over time.

## Should You Be Interested in the Tutorial?

### If You Are a University Lecturer
You may find this tutorial useful as either a standalone resource or a supplement to your own course material. The quickest way to assess its value is to examine the expected student outcomes after completing all exercises. These can be reviewed by:  
- Inspecting the **task boxes** at the beginning of each exercise (see the [latest release](https://github.com/Foly93/MD_FromBasicsToApplication/blob/main/releases/LiveCoMS_tutorial_v1.pdf)).  
- Reading **model reports** compiled by one of the authors during their master’s studies (available in the repository subdirectories, e.g., [01_INTRO](https://github.com/Foly93/MD_FromBasicsToApplication/tree/main/01_INTRO), [02_POLYMER_FF](https://github.com/Foly93/MD_FromBasicsToApplication/tree/main/02_POLYMER_FF), etc.).  
It’s also recommended to read the *Introduction* and *Scope* sections of the [latest release](https://github.com/Foly93/MD_FromBasicsToApplication/blob/main/releases/LiveCoMS_tutorial_v1.pdf).  


### If You Are Not a University Lecturer
We recommend starting with the first exercise and working through the course material step by step. Please note that the prerequisites listed in the [latest release](https://github.com/Foly93/MD_FromBasicsToApplication/blob/main/releases/LiveCoMS_tutorial_v1.pdf) must already be installed on your workstation, as neither the authors nor the tutorial provide installation guidance.  
Before starting, and especially if you are in doubt, we encourage you to briefly explain your situation and ask whether attempting the tutorial would be worthwhile for you. You can do this by opening an issue on GitHub or by sending an email to the corresponding authors.

## The Content in All Brevity
Molecular dynamics simulations are powerful because of their complexity: the method is both robust and versatile. At the same time, the multitude of choices a scientist must make can be overwhelming. Even experienced biophysicists sometimes struggle to follow best practices or to avoid common pitfalls. It is not unusual for professional researchers to design flawed simulations from time to time, which can lead to wasted time and effort when work must be redone.
This tutorial addresses these challenges by teaching the **foundations of MD simulations first**, then moving on to **applied methods and advanced analysis**, and finally guiding students toward **independent simulations** whose structure mirrors the workflow of computational biophysicists and physical chemists.  
**Foundations (Ex1–Ex5).** These exercises cover core concepts such as integration algorithms, timesteps, forcefield parameters, thermostats, barostats, and long-range electrostatics. Each topic is taught proactively: exercises are designed to yield real-world observables, which can then be compared against literature values to distinguish between good and poor algorithmic choices. This approach keeps learners engaged, encourages critical thinking, and trains them to test and evaluate their own setups - an essential skill for future projects. Special emphasis is placed on raising awareness of best practices and common mistakes to avoid.  
**Applied methods (Ex6–Ex7).** These exercises introduce more advanced calculations of free energy properties that connect molecular dynamics simulations with quantities of interest in experimental contexts. Students estimate solvation free energies and binding free energies at a simplified level, balancing time efficiency with conceptual clarity while also connecting to wet lab experiments.  
**Independent applications (Ex8–Ex9).** The final exercises are more self-directed, requiring students to rely on the knowledge and intuition gained from earlier work. First, they conduct and analyze a protein-in-solvent simulation based entirely on their own methodological choices. Finally, they perform a guided virtual screening exercise, which emphasizes managing and understanding the many software components involved in producing meaningful results.

### Versions

- [on github](https://github.com/Foly93/MD_FromBasicsToApplication/blob/main/releases/LiveCoMS_tutorial_v1.pdf)
- [on LiveCoMS](https://doi.org/10.33011/livecoms.6.1.3797)
     	 

## Get Involved

Help is always welcome, no matter the format. Any criticism, remarks, and participation are highly appreciated, and the most basic help is to work through the tutorials and ask questions that boggle your mind along the way.

### Submit a Pull Request

Suppose you are a lecturer with ideas on improving the exercises, tasks, or course material. In that case, proposing changes via a "pull request" is encouraged, allowing us to track your contributions. In case of substantial contribution, the author list in this article can be appended to subsequent versions of this work, as they would be for a software project. New versions of this work are assigned unique, cite-able DOIs and constitute preprints to be cited as interim research products.

## Authors
- Luis Vollmers (TUM Alumni)
- Shu-Yu Chen (ETHZ)
- Maria Reif (TUM Alumni)
- Tristan Alexander Mauck (TUM Alumni)
- Martin Zacharias (TUM)

Your name, too, can go here if you help us substantially revise/extend the paper.

### Citing This Work

```
@article{Vollmers2025,
author = {Vollmers, Luis and Chen, Shu-Yu, and Reif, Maria and Mauck, Tristan Alexander and Zacharias, Martin},
title = {Molecular Dynamics: From Basics to Application [Article v1.0]},
journal = {Living Journal of Computational Molecular Science},
volume = {6},
number = {1},
pages = {3797},
year = {2025},
doi = {11.33011/livecoms.6.1.3059},
}
```
