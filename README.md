# ATCM  
**A**lex **T**halamo-**C**ortical **M**odel — extended thalamo-cortical neural mass model for M/EEG / electrophysiology

---

## Table of Contents

1. [Overview](#overview)  
2. [Features](#features)  
3. [Model Structure & Theory](#model-structure--theory)  
4. [Getting Started](#getting-started)  
   * [Requirements](#requirements)  
   * [Installation](#installation)  
   * [Directory / File Structure](#directory--file-structure)  
5. [Using ATCM](#using-atcm)  
   * [Example Scripts](#example-scripts)  
   * [Running Simulations](#running-simulations)  
   * [Parameter Estimation / Inversion](#parameter-estimation--inversion)  
6. [Analyses & Outputs](#analyses--outputs)  
7. [Customization](#customization)  
8. [Contributing](#contributing)  
9. [Licence](#licence)  
10. [Citations](#citations)  

---

## Overview

ATCM is a MATLAB toolbox implementing an **extended thalamo-cortical neural mass model** (TCM) suitable for M/EEG and related electrophysiological data. It builds on canonical thalamo-cortical motifs, using Morris-Lecar–like dynamics and neural mass populations to model cortical laminar structure, thalamic relay/reticular interactions, and modulatory effects.  

It supports:

- time-series simulations under different inputs  
- spectral / transfer function analyses  
- modal decompositions and singular spectrum methods  
- parameter estimation (via e.g. variational Laplace within DCM frameworks)  

---

## Features

- Rich dynamic behavior (oscillations, spindles, gamma, etc.)  
- Flexible input types (deterministic pulses, continuous drives, noise etc.)  
- Multiple integration algorithms & time-stepping options  
- Tools for spectral / transfer function computation  
- Support for prior specification, modal decomposition, and fitting/inverting the model to data  

---

## Model Structure & Theory

- Populations included: cortical (e.g. superficial, deep), inhibitory interneurons, thalamic relay and reticular nuclei.  
- Neuronal dynamics via membrane potential / activation gating (Morris-Lecar-like) equations.  
- Synaptic coupling between populations, with both excitatory and inhibitory connections.  
- Optional modulation (e.g. neuromodulatory effects) in some versions / scripts.  
- Modal decomposition & spectral analyses rely on linearization / transfer-function approximations in certain regimes.  

---

## Getting Started

### Requirements

- MATLAB (tested on version 2023a)  
- Signal Processing Toolbox (for spectral analyses)  
- Optional: toolboxes for optimization / variational inference (e.g. fitVariationalLaplaceThermo.m in https://github.com/alexandershaw4/aLogLikeFit)  

### Installation

1. Clone this repository:  
   ```bash
   git clone https://github.com/alexandershaw4/atcm.git
