# R Calculator Project

This repository contains a menu-driven calculator written in R (`calc.R`). The script implements Mathematics, Descriptive Statistics, Numerical Methods, and Simulation tools accessible via interactive console menus.

## How to run

1. Open R or RStudio and set your working directory to the project folder, for example:

```r
setwd("/path/to/calculator_project")
source("calc.R")
# The calculator will start automatically and display the MAIN menu
```

The top-level function `MAIN()` launches the interactive menu. Choose a section and follow on-screen prompts.

## File

- `calc.R` — Main script with all functions and menus. The script calls `MAIN()` at the end to start the calculator.

## High-level sections and functions

The project is organized into four top-level sections, each exposed through the main menu. Below each section is a list of functions found in `calc.R`, a short description, input expectations, and example usage notes.

### 1) Mathematics

Functions:
- `AdditionSubtraction1()`
  - Purpose: Sum numbers entered interactively.
  - Inputs: Enter numbers (press ENTER after each). Finish by pressing ENTER twice.
  - Output: Prints `Result = <sum>` to console.

- `Multiplication2()`
  - Purpose: Multiply numbers entered interactively.
  - Inputs: Same pattern as addition.
  - Output: Prints `Result = <product>`.

- `Divide2Nos3()`
  - Purpose: Divide two numbers (a / b).
  - Inputs: Enter exactly two numbers.
  - Output: Prints `Result = <quotient>`. Validates division by zero.

- `Exponent4()`
  - Purpose: Compute base^power.
  - Inputs: Enter base then power.
  - Output: Prints `Result = <value>`. Handles 0^0 and division-by-zero for negative powers.

- `Factorial5()` and helper `fact(x)`
  - Purpose: Compute factorial of a non-negative integer.
  - Inputs: Enter a non-negative number (integer expected).
  - Output: Prints `Result = <factorial>`.

- `Permutation6()`
  - Purpose: Compute P(N, R) = N! / (N-R)!.
  - Inputs: Enter N and R (integers).
  - Output: Prints `Result = <permutation>`.

- `Combination7()`
  - Purpose: Compute C(N, R) = N! / ((N-R)! * R!).
  - Inputs: Enter N and R (integers).
  - Output: Prints `Result = <combination>`.

- `AP8()`, `GP9()`, `HP10()`
  - Purpose: Generate Arithmetic, Geometric, and Harmonic Progressions respectively.
  - Inputs: First term, common difference/ratio, and number of terms.
  - Output: Prints the sequence to console.

- `SimpleInterest11()` and `CompoundInterest12()`
  - Purpose: Compute simple and compound interest with formatted output.
  - Inputs: Principal, rate (% per annum), time (years), and compounding frequency for compound interest.
  - Output: Prints breakdown (Principal, Rate, Time, Interest, Total Amount).

- `QuadraticRoots13()`
  - Purpose: Solve quadratic ax^2 + bx + c = 0 and print real or complex roots.
  - Inputs: Coefficients a, b, c.
  - Output: Prints roots with proper formatting and validation for a != 0.

- `HCF14()` and helper `gr_cm_di(a, b)`
  - Purpose: Compute HCF (GCD) of two positive integers using Euclidean algorithm.
  - Inputs: Two integers.
  - Output: Prints `Result = <gcd>`.

- `LCM15()`
  - Purpose: Compute LCM of two positive integers.
  - Inputs: Two integers.
  - Output: Prints `Result = <lcm>`.

- Matrix operations: `MatrixAdd16()`, `MatrixSubtract17()`, `MatrixTranspose18()`, `MatrixMulti19()`, `MatrixInv20()`
  - Purpose: Basic matrix arithmetic: addition, subtraction, transpose, multiplication, inverse.
  - Inputs: Matrix dimensions and elements entered row-wise. Inverse requires a square, non-singular matrix.
  - Output: Prints resulting matrix.

- `Mathematics()`
  - Purpose: Interactive menu for all mathematics functions above. Choose items 1–21.

### 2) Descriptive Statistics

Functions:
- `Raw1()`
  - Purpose: Compute descriptive statistics from raw observations entered interactively.
  - Inputs: Enter observations (press ENTER after each, finish by pressing ENTER twice).
  - Output: Prints count, min, max, AM, GM (if possible), HM (if possible), quartiles, mode(s), range, IQR, mean deviations, variance, standard deviation, coefficient of variation, skewness/kurtosis measures.

- `Frequency2()`
  - Purpose: Compute descriptive statistics for frequency data.
  - Inputs: Enter distinct xi values (space-separated), then frequencies fi. Both input vectors must match length.
  - Output: Same statistics as `Raw1()` but computed for grouped data.

- `Bivariate3()`
  - Purpose: Compute bivariate summary (means, variances, covariance, correlation, regression Y on X, probable error and significance test).
  - Inputs: Enter Xi values then Yi values (press ENTER after each list).
  - Output: Prints bivariate statistics and regression equation.

- `DescStats()`
  - Purpose: Interactive menu for descriptive statistics (Raw Data, Frequency Data, Bivariate Data).

Notes: Input validation is present for non-numeric values and mismatched lengths.

### 3) Numerical Methods

Functions:
- `BisectionMethod1()`, `SecantMethod2()`, `RegulaFalsi3()`, `NewtonRaph4()`
  - Purpose: Root finding methods. Each prompts to choose one of three sample functions and then requests initial guesses/intervals.
  - Inputs: Choice of function and initial approximations.
  - Output: Iteration logs, final root, and a plot of the function with the root highlighted.

- `ForwardSub5()` and `BackwardSub6()`
  - Purpose: Solve triangular linear systems using forward/backward substitution. They accept the coefficient matrix and RHS vector.
  - Inputs: Order (n), n*n matrix entries (row-wise), and n RHS entries.
  - Output: Prints the solution vector.

- `Lagrange7()`
  - Purpose: Lagrange interpolation for given data points.
  - Inputs: Xi values, FXi values (matching lengths), and x to interpolate at.
  - Output: Prints interpolated value.

- `NewtonForward8()` and `NewtonBackward9()`
  - Purpose: Marked 'Under Construction'. No implementation currently.

- Numerical integration: `Trapezoidal10()`, `Simpson13rd11()`, `Simpson38th12()`
  - Purpose: Numerical integration using Trapezoidal, Simpson 1/3, and Simpson 3/8 rules. Prompts for function choice, limits a and b, and number of intervals.
  - Inputs: Choice of function, limits, number of intervals (n) with validity checks (e.g., Simpson 1/3 requires even n).
  - Output: Prints result and plots the function with shaded areas.

- `Numerical()`
  - Purpose: Interactive numerical methods menu.

### 4) Simulation

Functions:
- `BuffonNeedle1()`
  - Purpose: Approximate pi using Buffon's needle simulation.
  - Inputs: Number of iterations, needle length `l`, distance between lines `d`.
  - Output: Approximate value of pi.

- `MidSqRNG2()`
  - Purpose: Mid-square RNG using a 4-digit seed.
  - Inputs: Number of random numbers and a 4-digit seed.
  - Output: Prints sequence of generated integers.

- `CongRNG3()`
  - Purpose: Linear congruential RNG with checks for Knuth conditions where possible.
  - Inputs: Number of RNGs, multiplier `a`, increment `c`, modulus `m`, seed `x`.
  - Output: Prints sequence and whether Knuth conditions were met.

- `MultiRNG4()`
  - Purpose: Multiplicative congruential RNG (no increment).
  - Inputs: Number of RNGs, multiplier `a`, modulus `m`, seed `x`.
  - Output: Prints sequence.

- `Simulation()`
  - Purpose: Interactive menu for simulation methods.

## Notes and assumptions

- The script uses `scan()` and `readline()` heavily; it's intended for interactive console use (R terminal or RStudio console).
- Some functions print plots using base R plotting; these require an interactive graphics device.
- `NewtonForward8()` and `NewtonBackward9()` are placeholders and currently print "Under Construction".
- The script imports `DescTools` at the top; ensure that package is installed for certain helper functions (e.g., `Primes()` used in RNG checks).

## Example quick session

1. Start R in the project folder and run `source("calc.R")`.
2. From `MAIN` choose `1` for Mathematics.
3. Choose `1` to sum numbers; enter `1`, `2`, `3`, and press ENTER twice. The script prints `Result = 6`.

## Next steps / improvements

- Add unit tests for core functions (factorial, gcd, lcm, basic matrix ops).
- Implement `NewtonForward8()` and `NewtonBackward9()`.
- Replace `setwd(...)` at top with a relative path or remove; using `setwd` in scripts can be intrusive.
- Add argument-based wrappers so functions can be used programmatically (not only via interactive input).

---

Requirements coverage:
- Read `calc.R` and documented all functions: Done
- Created `README.md` with function names and usage: Done

Generated on 2025-09-19.
