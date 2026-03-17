# Satellite Orbital Dynamics Modelization

Development Period: Apr 2024 - Jun 2024 (Original C), Mar 2026 (C++ Refactoring)
Languages: C, C++ (C++17)

## Description

This repository contains a numerical solver that models the orbital motion of a satellite in a central gravitational field. The mathematical core integrates the equations of motion using a fourth-order Runge–Kutta (RK4) method, allowing the analysis of orbital evolution, energy conservation, and variable-mass dissipative effects (powered flight dynamics).

The project is split into two distinct implementations:

* **C Version (`satellite.c`)**: The original implementation. It is functionally and mathematically correct, utilizing standard C procedural paradigms.
* **C++ Version (`satellite.cpp`)**: A modernized C++17 refactoring optimized for performance, safety, and maintainability. Key features include:
    * **Zero-Allocation Fast-Path**: The mathematical integration loop is strictly decoupled from disk I/O and heap allocations to guarantee stable execution times.
    * **Contiguous Memory Management**: Uses pre-allocated `std::vector` capacity to prevent heap fragmentation and maximize CPU cache hits.
    * **Safety & RAII**: File descriptors and memory boundaries are automatically managed via RAII. Control flow is enforced via strict typing (`enum class`).

## How to use

Both versions are designed as standalone files for immediate compilation.

**To compile the C++ version (Recommended):**
```bash
g++ -O3 -std=c++17 satellite.cpp -o satellite_cpp.x
./satellite_cpp.x
```
**To  compile the C version:**
```bash
gcc satellite.c -o satellite.x -lm
./satellite.x
```
