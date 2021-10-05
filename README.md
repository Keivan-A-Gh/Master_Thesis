# Release of High-Pressure (70 MPa) Hydrogen into the Air
**By Keivan A.Gh.**
<p align="center">
<img width="800" height="222" src="https://user-images.githubusercontent.com/44712667/135481464-23cf1baf-371a-4cc0-b585-5b83efad485f.gif">
</p>

## Abstract of This Project
<p>Hydrogen as an energy carrier is one of the available alternatives to fossil fuels for decarbonizing the global energy system. Regarding to the very low volumetric energy content of gaseous hydrogen, for practical and economical storage, it has to be either pressurized at ambient temperature (Compressed Gaseous Hydrogen – CGH2) or liquefied at cryogenic temperatures (Liquid Hydrogen – LH2). On the subject of CGH2, using high-pressure hydrogen reservoirs from 25 MPa up to 70 MPa is conventional; and among these, utilizing 70-MPa hydrogen reservoirs as fuel cells in transportation or as storage tanks, for instance, in local hydrogen refueling stations is the common approach.</p>
<p>Having a leakage in such a high-pressure reservoir will form shock waves in front of the released hydrogen causing a temperature rise which may be intensified by existing obstacles in that environment; the presence of these obstacles and confinement may also enhance the hydrogen-air mixing. Because of the wide flammability range of hydrogen, this might lead to its spontaneous ignition if these high-temperature, well-mixed regions could last as long as the hydrogen-air induction time. In this thesis, it is tried to numerically investigate the possibility of this kind of scenario.</p>
<p>In order to simulate a release of 70-MPa hydrogen into the air (treated as a dual pseudo species), in the initial attempts, the USN-FLIC code was tried, but because the results were not convincing, it was decided to use the OpenFOAM software as an alternative. Considering the restrictions of the solvers of OpenFOAM (v.7), combinations of solvers along with new thermophysical models are used to be able to overcome the so-called Riemann problem of shock waves in a non-ideal, multi-component, non-reacting mixture. Furthermore, to validate this method, the shock-tube problem is solved and the results are compared with available data of similar cases.</p>
<p>According to the results, potentially hazardous regions are formed in the domain that are mostly related to the interaction of the reflected, rarefaction, and normal shock waves inside the flow field. Although there are some inconsistencies between the results of the simulations in estimating the flow properties, generally, the risk of hydrogen auto-ignition in these regions is quite high. But for having a better understanding about the effects of distances between wall boundaries, it is shown that the simulation should be done in three dimensions.</p>

## Introduction to MASTER_THESIS Repository
This repository contains 2 main folders:

  1. **[Documents](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/Documents):**
      - You can find the [report](./Documents/Master's%20Thesis%20-%20Final%20Edition.pdf) of this project and some [visualized results](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/Documents/Visualization%20of%20the%20Results) in this folder.
  2. **[OpenFOAM v.7](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7):**
      - This folder contains all the implemented [OpenFOAM utilities](./OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models) as well as the new generated [thermophysical coefficients](./OpenFOAM%20v.7/New%20Thermophysical%20Coefficients) and the [initial setup](./OpenFOAM%20v.7/Tutorial%20(OF%20Simulation%20Setup)) of the OpenFOAM simulation.

***N.B.***  For further explanations about each section, you can refer to the README files of these folders.


## Citation
If you want to refer to this project or its corresponding implemented utilities, please cite the following work:
  - [![DOI](https://zenodo.org/badge/411550302.svg)](https://zenodo.org/badge/latestdoi/411550302)
