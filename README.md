# Numerical Analysis Methods Implementation

This repository contains implementations of various numerical methods learned during my Numerical Analysis course. It serves as both a practical demonstration of the concepts and a reference for future use.

## Course Topics & Implementations

### 1. Linear System Solvers
- **Conjugate Gradient Method** (`CG.m`)
  - Iterative method for solving sparse linear systems
  - Particularly efficient for symmetric, positive-definite matrices

- **Gauss-Seidel Method**
  - Iterative solution for linear equation systems
  - Improved convergence over Jacobi method

- **Successive Over-Relaxation** (`q2_SOR.m`)
  - Enhanced version of Gauss-Seidel method
  - Introduces relaxation parameter for faster convergence

### 2. Root Finding Methods
- **Fixed-Point Iteration** (`q1.m`)
  - Simple iterative method for finding function roots
  - Implementation includes convergence criteria

- **Newton's Method**
  - Quadratic convergence for well-behaved functions
  - Requires function derivative

### 3. Eigenvalue Computation
- Multiple implementations exploring different approaches:
  - Built-in `eig` function usage
  - QR decomposition method
  - Applied to test matrices A, B, and C

### 4. Differential Equations
- **Runge-Kutta Method** (`Runge-Kutta.m`)
  - Fourth-order implementation
  - Handles initial value problems
  - Higher accuracy compared to simpler methods

### 5. Numerical Integration
- **Trapezoidal Rule** (`Trapezoid.m`)
  - Implementation for definite integral approximation
  - Error analysis included

## Project Structure
- Core implementation files (`.m` files)
- Test cases and example problems
- Supporting equation systems (`equation.m`)
- Detailed documentation (`.docx`)

## Learning Outcomes
Through this project, I've gained practical experience in:
- Implementing various numerical methods
- Understanding convergence properties
- Error analysis and method selection
- MATLAB programming for mathematical applications

## Usage
Each implementation includes:
- Function documentation
- Example usage
- Error handling
- Performance considerations

For detailed explanations and mathematical background, refer to the accompanying documentation.
