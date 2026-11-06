WOLBACHIA ABM — INPUT DATA
-------------------------
This README describes the environmental input files and calibrated burn-in snapshot used in the mosquito agent-based model.

1. ambiental.csv — SIMULATION TEMPERATURE DATA
----------------------------------------------
Daily temperature records for Niterói (2017–2020), from which Zone 3 temperatures are used in the mosquito population simulations.

2. ambiental16.csv — BURN-IN TEMPERATURE DATA
--------------------------------------------
Daily temperature records for Niterói (2016), from which Zone 3 temperatures are used during the 364-day population burn-in period.

3. K25000.jls — CALIBRATED BURN-IN SNAPSHOT
------------------------------------------
Serialized mosquito population state generated after burn-in using the calibrated carrying-capacity parameter K = 25,000, selected through simulated minimum-distance calibration against observed Wolbachia prevalence.

The saved burn-in state is used to initialize the subsequent simulation experiments.




