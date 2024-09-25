Remarks/Observations: 
- Coarse mesh was found to converge to a solution without the need of the cw1 term ramping function, however near-wall results (e.g. velocity and eddy viscosity) differs "somewhat" from the literature found on https://turbmodels.larc.nasa.gov/flatplate_sa.html. 
- Fine mesh was found to not converge (i.e. resulting in NANs residuals) regardless of the use of the cw1 term ramping function, however with the ramping function it does run for longer.
- The log file and results (in the "results" folder) were obtained with the coarse mesh.
