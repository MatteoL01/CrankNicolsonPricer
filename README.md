# Option Pricing and Financial Models

This repository implements financial option pricing models, including the Crank-Nicolson finite difference method, and the Black-Scholes model. It includes various utilities and algorithms to calculate option prices, Greeks, and interest rate interpolations.

---

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)

---

## Overview

This project provides tools for:
- Pricing European and American options using the Crank-Nicolson method.
- Computing the Greeks (\(\Delta\), \(\Gamma\), \(\Theta\), \(\rho\), \(\nu\)).
- Using the Black-Scholes model as a benchmark.
- Handling piecewise interest rate structures with linear interpolation.

It supports both user-defined parameters and default settings, making it adaptable for financial modeling tasks.

---

## Features

- **Crank-Nicolson Method**: Solves the partial differential equation governing option pricing.
- **Black-Scholes Model**: Calculates exact option prices for European options.
- **Utility Functions**:
  - Cumulative and probability density functions for the normal distribution.
  - Linear interpolation of interest rates.
  - Boundary condition handling for numerical methods.
- **Output**: Results are exported to a text file for easy analysis.
- **Modular Design**: Clear separation of functionality with reusable components.

---


### Installation
Clone the repository:
   ''bash
   git clone https://github.com/yourusername/option-pricing.git
   cd option-pricing

### Usage

g++ -std=c++17 -o option_pricing main.cpp Matrix.cpp Diffmethod.cpp utils.cpp

