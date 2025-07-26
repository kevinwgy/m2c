# Spherical explosion simulation (Stage 2)

The file `SphericalShock.txt` contains pressure, velocity, and 
density profiles generated from simulating Stage 1 (material 
ignition and detonation).

The fluid state within the spherical explosive region is 
initialized by interpolating values from this file. The file is
provided to M2C through the `UserDataFile` keyword under the 
initial conditions.

# Running the example case

To run the simulation, use:

```sh
../../../m2c input.st 2>&1 | tee log.out
```

If you have access to multiple CPUs, you can perform the simulation
in parallel. For example, to it on eight processors, use:

```sh
mpiexec -n 8 ../../../m2c input.st 2>&1 | tee log.out
```

For execution using a job scheduler such as SLURM, you can use the 
example bash script run.sh, which can be submitted with:

```sh
sbatch run.sh
```

Before launching the simulation, ensure that the bash script is 
properly configured with the correct user ID, group ID, and requested 
computational resources.
