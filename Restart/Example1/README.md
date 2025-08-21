# Introduction

This example demonstrates how to restart a simulation using the solution from a previous run as the initial condition. The physical setup is identical to `Test1` in `Tests/RiemannProblems`. It is a one-dimensional two-material Riemann problem with two constant initial states for x<0.5 and x>0.5, respectively. As time advances, the solution features a shock wave moving to the right and a rarefaction fan moving and expanding to the left.


### Directory Structure

This example contains three subdirectories:

- **`InitialRun`**  
  Runs the original simulation (`Test1`), but with a shorter final time (`0.1` instead of `0.2`).
  
- **`Restart`**  
  Restarts the simulation from the final state of `InitialRun` and continues for another `0.1` time units.
  
- **`RestartMod`**  
  Similar to `Restart`, but starts from a **modified version** of the `InitialRun` final solution.



# Instructions

### 1. Run the initial simulation

```bash
cd InitialRun
mpiexec -n 2 [path-to-m2c-executable] input.st
```

The solution files will be written to the `results` directory.

---

### 2. Prepare the restart files

Navigate to the `Restart` directory and copy the final solution from the initial run:

```bash
cd ../Restart
cp ../InitialRun/results/solution_0010.vtr IC
```

Then, compile the user-defined initial condition code in the `IC` folder. Edit `CMakeLists.txt` to provide the path to the M2C source code.
> **Note:** This requires VTK version 9.3 or later.

```bash
cd IC
cmake .
make
```

---

### 3. Run the restart simulation

From the `Restart` directory, run:

```bash
mpiexec -n 2 [path-to-m2c-executable] input.st
```

The `input.st` file is nearly identical to that of the initial run, except it specifies the **user-defined initial condition**.

Notice that the input file `input.st` is almost identical to that in the initial run, except for the specification of the user-defined initial condition.

---

### 4. Run the modified restart simulation

Enter the `RestartMod` subdirectory. Repeat **Step 2** above. Then, go to the `USDL` (**USer-Defined Solution**) subdirectory and compile the solution modification file:

```bash
g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedSolution.h -c UserDefinedSolution.cpp; g++ -shared UserDefinedSolution.o -o UserDefinedSolution.so; rm UserDefinedSolution.o
```

Then, return to `RestartMod` and run:

```bash
mpiexec -n 2 [path-to-m2c-executable] input.st
```

> **Note:** The input file here is almost identical to that in `Restart`, except for the specification of **"PrescribedSolution"**.

---
