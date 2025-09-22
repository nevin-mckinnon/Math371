# Optimizing Cancer Radiation Therapy via Mathematical Modeling
This repository contains the MATLAB code for a final project for Math 371. The project's goal is to develop a mathematical model of tumor growth and its interaction with healthy tissue, and then use this model to find an optimal radiation therapy plan. The optimization aims to maximize the destruction of tumor cells while minimizing damage to surrounding healthy cells.

## Project Overview
Radiotherapy is a cornerstone of cancer treatment, but it presents a critical challenge: how to administer a dose strong and long enough to eradicate a tumor without causing unacceptable harm to the patient. This project addresses that challenge by creating a computational simulation that can test and optimize treatment strategies virtually.

The model is built on a system of Ordinary Differential Equations (ODEs) and is simulated in three distinct phases:

Pre-Treatment: Natural growth and competition between cell populations.

During Treatment: The effects of radiation on both tumor and healthy cells.

Post-Treatment: The long-term outcome, observing recovery and potential tumor regrowth.

The core of the project is an optimization algorithm that searches for the ideal radiation dose (D) and duration (trad) to achieve the best clinical outcome within the simulation.

## The Mathematical Model
The simulation is primarily based on two well-established biological and radiobiological models:

### 1. Cell Growth and Competition (Lotka-Volterra Model)
Instead of simple exponential growth, the model implements logistic growth, which assumes that a cell population's growth slows as it approaches a carrying capacity (K). Furthermore, it incorporates competitive interactions between the tumor (VT) and healthy (VH) cell populations, based on the Lotka-Volterra competition model. This creates a more realistic dynamic where both cell types compete for shared resources.

### 2. Radiation Effects (Linear-Quadratic Model)
The model simulates cell death from radiation using the Linear-Quadratic (LQ) model. This is a standard radiobiological model that relates the fraction of cells surviving a dose of radiation (D) to two key parameters: alpha (representing single-hit, irreparable damage) and beta (representing two-hit, repairable damage). This allows the simulation to accurately model how both tumor and healthy cells respond differently to radiation.

## Methodology & Algorithms
The project executes a sophisticated pipeline to determine and simulate the optimal therapy.

### Simulation Phase 1: Before Radiation
The simulation begins by modeling the initial state of the system. It solves the system of ODEs to show how the tumor and healthy cell populations evolve and compete before any medical intervention.

### Optimization of Treatment Parameters
This is the heart of the project. Using the state of the system right before treatment as a starting point, the script employs a powerful optimization algorithm to find the best treatment parameters.

**Algorithm:** The optimization is performed using MATLAB's fmincon function, configured to use the Sequential Quadratic Programming (SQP) algorithm.

**Objective:** The goal is to find the dose (D) and duration (trad) that minimizes an objective function, which is designed to reward the reduction in tumor volume while penalizing the reduction in healthy cell volume.

### Simulation Phase 2 & 3: During and After Radiation
With the optimal parameters identified, the script runs the simulation through the treatment period and into a post-treatment phase. This demonstrates the full effect of the optimized plan, showing the decline of the tumor, the temporary damage and subsequent recovery of healthy tissue, and the final state of the system.

For solving the differential equations in all simulation phases, the model uses MATLAB's ode23 solver, which is based on the efficient Bogacki-Shampine method.

# How to Run the Code
  1. Ensure you have MATLAB installed.
  2. Clone this repository and open the main .m script file.
  3. The parameters for the model (growth rates, competition factors, etc.) can be adjusted in the initial "Constants and Initial Conditions" section.
  4. Run the script.

# Results
The script will produce a stacked area plot visualizing the volume of tumor cells, non-dividing (dead) cells, and healthy cells over the entire simulation time. The plot is annotated with vertical lines indicating the start and end of the radiation treatment, and a text box displays the calculated optimal dosage and treatment duration.
