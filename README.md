# TPC analysis code

## Introduction

This repositroy is a collection of code written for the project "First Conceptual Design and Studies for a Tracking Time Projection Chamber for the Belle II Experiment", which is a Master's thesis written by Andreas Löschcke Centeno at the University of Bonn in the group of Jochen Dingfelder under supervision of Peter Lewis and Christian Wessel. 
The thesis can be found [here](https://docs.belle2.org/record/2631/files/BELLE2-MTHESIS-2021-073.pdf).

## Project

This project was entirely simulation based and the main code for the simulation can be found in the Belle II software repository. 

The files here are a loose collection of analysis code written for the TPC, which was not appropriate to be added to the Belle II repository. It is saved here in case someone else wants to pick up the TPC project. The files here mainly consist of jupyter notebooks used for visualizing and analysing properties of the TPC. 

Additionally, there is one file (`eventoverlay.cc`) which is a ROOT macro used to create the event overlap in a TPC. In the folder 'manual' there is a short guide to this overlay code. The guide is intended to be read alongside the code.



## Contact

Andreas Löschcke Centeno:
- private: a.loeschcke@gmail.com 
- work: a.loeschcke-centeno@sussex.ac.uk

Peter Lewis:
- work: lewis@physik.uni-bonn.de

Christian Wessel:
- work: wessel@physik.uni-bonn.de