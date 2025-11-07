# Variational Monte Carlo for the Quantum Harmonic Oscillator (QHO)

This project is from independent research and development conducted during my **Master of Science program in Computational Physics**. The entire methodology, analysis, and implementation of all numerical solvers were executed solely by me, and the full work is documented in the accompanying report and code.

---

## Table of Contents

* [About The Project](#about-the-project)
* [Core Objectives](#core-objectives)
* [Languages and Libraries](#languages-and-libraries)
* [Methods Implemented](#methods-implemented)
* [Key Findings](#key-findings)
* [Getting Started](#getting-started)
* [Full Project Report](#full-project-report)
* [Contact](#contact)

---

## About The Project

This work focuses on using the **Variational Monte Carlo (VMC)** method coupled with the **Metropolis algorithm** to determine the ground state energy of a system of one and two interacting electrons confined in a two-dimensional harmonic oscillator potential.

The project demonstrates how to use a variational parameter within a trial wave function (incorporating a Jastrow factor) to accurately model and capture the complex **electron-electron correlation effects** in the system.


## Core Objectives

1.  **Implement VMC:** Develop a robust Variational Monte Carlo simulation using the Metropolis algorithm.
2.  **Optimize Trial Function:** Employ a **Golden Section Search** to efficiently minimize the system's energy variance by finding the optimal variational parameter ($\alpha$).
3.  **Verify One-Electron System:** Accurately compute the ground state energy of a single electron in a 1D harmonic oscillator for validation against known analytical solutions.
4.  **Analyze Two-Electron Correlation:** Apply the VMC method to a two-electron system with varying Coulomb repulsion strength ($\lambda$) to quantify the effects of electron-electron interaction.

---

## Languages and Libraries

| Category | Tools & Libraries | Competency Demonstrated |
| :--- | :--- | :--- |
| **Language** | Python | Efficient development and handling of statistical and numerical algorithms. |
| **Numerical** | NumPy, SciPy | Advanced array manipulation and implementation of the Metropolis algorithm and Golden Section Search. |
| **Visualization** | Matplotlib | Generating high-quality physics focused data plots. |

---

## Methods Implemented

The computational core of the project implements the following techniques:

| Method | Role in Project | Key Implementation Detail |
| :--- | :--- | :--- |
| **Metropolis Monte Carlo** | Core simulation engine. | Generates the ensemble of electron positions based on the probability distribution $|\Psi_T|^2$ to compute the expectation value of the local energy ($E_L$). |
| **Variational Principle** | Theoretical foundation. | Uses a trial wave function ($\Psi_T$) with a variational parameter ($\alpha$) to approximate the ground state energy by minimizing the variance. |
| **Jastrow Factor** | Modeling electron correlation. | Incorporated into the two-electron trial wave function to accurately account for the singular nature of the Coulomb interaction, particularly at strong repulsion ($\lambda$). |
| **Golden Section Search** | Optimization algorithm. | An efficient method used to iteratively find the optimal value of the variational parameter ($\alpha$) that yields the minimum variance/energy. |

## Key Findings

* **Validation Success:** The VMC successfully computed the ground state energy for the non-interacting one-electron QHO system, matching the analytical result $E=3\hbar\omega=3$ with high precision.
* **Correlation Capture:** For the two-electron system, the VMC accurately computed the expected energy $E=2$ when $\lambda=0$ (non-interacting). As the Coulomb repulsion ($\lambda$) increased, a clear minimum in the energy/variance plots was found, confirming the trial wave function's ability to capture complex electron-electron correlation (e.g., optimal $\alpha \approx 0.647$ for $\lambda=8.0$).
* **Variance as Metric:** The minimum variance was shown to reliably indicate the optimal variational parameter, reinforcing the theoretical strength of VMC.

---

## Getting Started

### Execution

To run the simulation and generate the results and visualizations, execute the core solver script:

```bash
python VMC_QHO.py
```
---

## Full Project Report

For a complete breakdown of the theoretical derivations and full results, please see the final project report:

[**Full Project Report (PDF)**](VariationalMC_QHO.pdf)

---

## Contact

I'm happy to hear your feedback or answer any questions about this project!

**Author** Rama Khalil 

**Email**  rama.khalil.990@gmail.com
