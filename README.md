# C++ Monte Carlo Simulation of 2D Ising Model
Monte Carlo simulation of 2D Ising Model on 1600+ spin node lattice in C++. Aim is to model and describe phase transition of magnetic substance through Monte Carlo sim of spin lattice, and determine how size of lattice influences behavior. Emphasis on optimized code and good memory/cache handling, escpecially in core loop (MonteCarloSim::NewState).

Simulation was conducted on 2d lattices of sizes 10x10, 40x40 and 100x100, achieving equilibrium after at most 5000 sweeps (timesteps) for most volitile temperature parameter. Longer time to settle for larger lattices.

Organized through IsingSubstance class, handling lattice and physical properties, and MonteCarloSim class handling the simulation itself.

Project handles measuring ergodic properties, such as heat capacity, as a function of temperature through TemperatureRun class.

Introduced child class for Ising Materials with impurities (fixed magnetic moment/spin), and solved energy optimum for these by simulated annealing.
