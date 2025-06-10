# PEM Fuel Cell Stack Simulation with Altitude Effects

This repository contains a Python simulation of a **Proton Exchange Membrane (PEM) Fuel Cell Stack**, which models how altitude affects the performance of a fuel cell. The model incorporates variations in **temperature**, **oxygen partial pressure**, and **voltage losses** due to changes in **altitude**. The simulation is inspired by the IMEKO 2007 paper:  
**"Analysis of Individual PEM Fuel Cell Operating Parameters for Design of Optimal Measurement and Control Instrumentation."** (Ref: Živko, Davor, and Vedran Bilas. "Analysis of individual PEM fuel  cell operating parameters for design of optimal measurement and control instrumentation." In 15th IMEKO TC4 Symposium on Novelties in Electrical Measurement and Instrumentation. 2007.).

---

## Features

- Models **temperature decrease** with altitude using the standard atmosphere lapse rate
- Adjusts **oxygen partial pressure** based on altitude
- Calculates the **reversible voltage** from thermodynamic principles
- Includes **activation**, **ohmic**, and **concentration** voltage losses
- Plots **Fuel Cell Stack Voltage vs. Current** at multiple altitudes
- Avoids common numerical issues (e.g., `log(0)`, division by zero)

---

## Output

The output is a graph showing the **voltage-current characteristics** of a PEM fuel cell stack at different altitudes (Limited to troposphere)

---

## Assumptions
- Parameters like 𝑖, 0, 𝛽, 𝑘, membrane thickness, diffusion layer conductivity, etc., are simplified based on typical values mentioned in paper. 
- Hydrogen partial pressure is assumed equal to ambient pressure, adjustments could be done for the pressurized H2 supply
- Water vapor partial pressure is estimated for 50% RH.

---

## Requirements

This project uses standard Python libraries:

- `numpy`
- `matplotlib`


