# jIIR2HW

This repository contains scripts and julia packages that permit to design second-order IIR filters for multiplierless hardware.

**Context**: multiplierless implementations of products by constants using shift-and-add architectures.

**Problem**: determine the coefficients of an IIR filter in fixed-point arithmetic such that the total cost of the shift-and-add architecture is minimal in terms of number of adders.

**Techniques**: our approach is based on a formalization of the problem as an instance of a Mixed Integer Linear Programming problem, which is subsequently solved using a third-party solver.

**Input**:
* frequency specifications
* wordlength

**Output**:
* filter coefficients
* adder graph structure for the optimal solution

**Guarantees** (in case of successful exit):
* optimality of the proposed solution
* a posteriori validation of frequency specifications

**Options**:
* bounded/unbounded adder depth
* avoid internal shifts in the adder graph structure
