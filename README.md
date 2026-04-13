# C++ Monte Carlo Simulation of 2D Ising Model
Monte Carlo simulation of 2D Ising Model on 1600+ spin node lattice in C++. Aim is to model and describe phase transition of magnetic substance through Monte Carlo sim of spin lattice, and determine how size of lattice influences behavior. Emphasis on optimized code and good memory/cache handling, escpecially in core loop (MonteCarloSim::NewState).

Ising materials show radically different time developments at a function of temperature, with signal volitility spiking at the critical temperature T_c approx 2.27 K - where variance goes to infinity. At this temperature, simulations take long to reach equilibrium (memoryloss of initial conditions), requiring highly efficient simulations. State of the art approach for the volitile state is typically the Wolff algorithm, where clusters of spins are flipped instead of single spins. The latter approach, single spin at a time, is the one implemented here, through a standard Metropolis algorithm using probabilities derived from statistical physics. The lack of algorithmic efficiency in our approach requires the code to be highly optimized in order to deal with temperatures close to T_c. 

Simulation was conducted on 2d lattices of sizes 10x10, 40x40 and 100x100, achieving equilibrium after at most 5000 sweeps (timesteps) for most volitile temperature parameter, taking longer time to settle for larger lattices. The project was a success: the scholastic yardstick of obtaining efficient timeevolution for 40x40 lattices was cleared, and the codebase manages to get interesting results for 100x100 and beyond efficiently.

Code relating to the lattice itself and its physical properties were organized through an IsingSubstance class. A MonteCarloSim class handled the simulation. In order to observe physical quantities (mainly heat capacity) over different temperatures, a TemperatureRun class was implemented. The project includes a child class for IsingSubstance, ImpureIsingSubstance, for simulating impurities in the material (fixed magnetic moment/spin, which cannot flip). Such impurities create complex boundary conditions for energy optima - the class SimulatedAnnealing uses exactly simulated annealing to find the optimized configurations across the lattice state space. The nice visuals of such optima are shown at the bottom of the readme.



Time evolution of energy (blue) and magnetization (red) for 100x100 lattice over 15000 timesteps, at T=2.3 K(close to T_c). Runtime: 2.59428s, using -03. We see convergence to equilibrium after <4000 timesteps.
<img width="560" height="420" alt="timeseries100_15k" src="https://github.com/user-attachments/assets/90f1fb7f-64ad-4b3f-bd2b-9a6c1fad4f5b" />

Time evolution of mean normalized energy(blue) and magnetization(red) for 40x40 lattice over 200000 timesteps, T=2.3 K. On larger timescales, the chaos close to critical temperature becomes appearent.


<img width="560" height="387" alt="time series" src="https://github.com/user-attachments/assets/f1be6b3d-fc36-4160-9b36-cb7a02b63386" />


Heat capacity as a function of temperature for lattice sizes 100x100 (blue), 40x40 (red), 10x10 (green). We discover that heat capacity scales with size of lattice.

<img width="560" height="420" alt="heatCap_102040" src="https://github.com/user-attachments/assets/0e235a06-2abe-42af-8133-cab326435134" />




We distribute impurities (spins/magnetic poles which cannot flip) randomly throughout the lattice with different concentrations, creating complex energy minema. Using simulatied annealing, starting at T=75 and cooling close to 0, we obtain optimized configurations! As expected, the configurations show larger "blocks" of the lattice pointing same direction (saving energy), broken up by fixed spins. The energy optimum of a pure substance would have all the spins pointing in the same direction, either up or down.

Energy optimum for 40x40 lattice with 0.09 share randomly distributed impurities:
<img width="560" height="420" alt="simulated annealing 40x40_2" src="https://github.com/user-attachments/assets/c3c81521-1362-4ec2-ab24-a63bea97c5c4" />

Energy optimum for 40x40 lattice with 0.40 share randomly distributed impurities:
<img width="560" height="420" alt="sim_annealing_40x40_0 4" src="https://github.com/user-attachments/assets/378027a4-2dc8-40a2-be30-a41674d56b37" />

Code for the different parts of the project has been /**/-ed out in main(). Feel free to try the code, although the graphics library is somewhat elaborate, and the better way might be to output a csv and visualize it somewhere else.
