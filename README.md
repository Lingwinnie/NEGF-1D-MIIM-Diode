# 1D NEGF Simulation of MIIM Diodes for Optical Rectennas

## Overview
This repository contains a MATLAB/Octave script developed to simulate the quantum transport properties of Metal-Insulator-Insulator-Metal (MIIM) diodes. It was originally built to evaluate the electrical feasibility of optical rectennas for thermophotovoltaic (TPV) applications. Following the geometric optimization of the multilayer structures for maximum light absorption, this code verifies the electrical constraints by calculating the tunneling current and the rectification ratio.

## Physical Model
The simulation relies on the Non-Equilibrium Green's Function (NEGF) formalism to model 1D ballistic quantum transport. The computational sequence includes:
* **Band Diagram Construction:** Establishes the equilibrium potential profile using the work functions of the metal electrodes and the electron affinities of the insulating layers.
* **Hamiltonian Formulation:** Uses a tight-binding approximation with finite differences to construct the device Hamiltonian and computes the self-energy matrices for the semi-infinite open contacts.
* **Equilibrium Transport:** Inverts the matrix to extract the retarded Green's function. It calculates the transmission probability via the Landauer formula and evaluates the spatially resolved Local Density of States (LDOS).
* **Current-Voltage (I-V) Characteristics:** Executes a voltage sweep. The potential profile is dynamically tilted for each bias point. The total current density is computed by integrating the transmission weighted by the difference in Fermi-Dirac distributions of the source and drain contacts.

## Physical Limitations
As a first-order engineering sizing tool, this model relies on explicit physical simplifications:
* **Constant Effective Mass:** The calculation assumes a constant average effective mass across the entire domain. Discontinuities and variations in effective mass at the heterojunctions are ignored.
* **Static Potential Profile:** The voltage drop under bias is approximated using an empirical capacitive divider model based on the relative permittivities of the respective oxides. It does not perform a self-consistent resolution with the Poisson equation.
* **1D Ballistic Approximation:** The formulation is strictly one-dimensional. It neglects transverse modes and all inelastic scattering mechanisms, such as electron-phonon or electron-defect interactions.

## Usage
The source code consists of a single script designed for direct execution in MATLAB or GNU Octave. Physical and geometric parameters including oxide thicknesses, potential barriers, relative permittivities, and voltage sweep ranges are defined at the beginning of the file and must be adjusted manually to simulate different material configurations.

## Author

**[Lingwinnie]** Master's Student in [Nanosciences and Nanotechnologies: Nanoscale and Quantum Engineering] at [Aix-Marseille Université].

