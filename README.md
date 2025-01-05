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

This project provides a comprehensive framework for option pricing and financial modeling, utilizing an Excel-based user interface for streamlined input handling and output visualization. The key features include:

- **Excel-Based User Interface**  
  Users can input option parameters (e.g., strike price, volatility, maturity) directly into an Excel file. The calculated results, including option prices and Greeks (\(\Delta\), \(\Gamma\), \(\Theta\), \(\rho\), \(\nu\)), are displayed in Excel for easy analysis and interpretation.

- **Crank-Nicolson Method**  
  Implements the finite difference method to accurately price European and American options.

- **Black-Scholes Model**  
  Provides a benchmark by calculating exact prices for European options using the Black-Scholes formula.

- **Piecewise Interest Rate Handling**  
  Supports scenarios where interest rates change over time by using linear interpolation to model the rate structure.

---

## Features

- **Crank-Nicolson Method**: Solves the partial differential equation governing option pricing.
- **Black-Scholes Model**: Calculates exact option prices for European options.
- **Utility Functions**:
  - Cumulative and probability density functions for the normal distribution.
  - Linear interpolation of interest rates.
  - Crout algorithm to solve the matrix system.
- **Output**: Results are exported to an Excel sheet.

---


## Installation
Clone the repository:
   ''bash
   git clone https://github.com/yourusername/option-pricing.git
   cd option-pricing

## Usage
Run the following code on the terminal:
    g++ -std=c++17 -o option_pricing main.cpp Matrix.cpp Diffmethod.cpp utils.cpp

## Project Structure
'''option-pricing/
│
├── include/
│   ├── Option.h            # Header for the general option class
│   ├── Diffmethod.h        # Header for finite difference methods
│   ├── Matrix.h            # Matrix operations header
│   └── utils.h             # Utility functions
│
├── src/
│   ├── Option.h            # Implementation of the general option class
│   ├── Diffmethod.cpp      # Implementation of finite difference methods
│   ├── Matrix.cpp          # Implementation of matrix operations
│   └── utils.cpp           # Implementation of utility functions
│
├── main.cpp                # Main entry point for the program
└── README.md               # Project documentation


