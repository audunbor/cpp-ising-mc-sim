#pragma once
#include <Eigen/Dense>
#include "IsingSubstance.h"
#include <random>
#include <iostream>
#include <matplot/matplot.h>
#include <memory>


namespace IsingProject{
    class MonteCarloSim {
        int n; //Number of sweeps
        int time;
        void newState();
        void timeStep();
        float autoCorrolation();
        int lag = 200;
        std::random_device rd; 
        std::mt19937 gen{rd()}; 
        std::uniform_real_distribution<float> dis{0.0f, 1.0f};
        std::uniform_int_distribution<> distr; 
        //Eigen::ArrayXXi findIndices() {return 2 - is.lattice*is.sumOfNeighbors/2;}
        //Eigen::ArrayXXf findProbabilities () {return findIndices().unaryExpr([this](auto i){return is.pRatios[i];});} //Find probability ratio of spins in lattice flipping 
        public:
        const Eigen::Array<float, 5, 1> pRatios; //Ratios of probabilities in Metropolis algorithm for state with {4,2,0-2,-4} sum of neighbor spins
        std::shared_ptr<IsingSubstance> is;
        MonteCarloSim(float T, int L, int n) : n{n}, time{0}, is{new IsingSubstance(T, L)}, distr{0, L-1}, pRatios{
            (Eigen::Array<float, 5, 1>{8,4,0,-4,-8}/(-T)).exp()
        } {};
        MonteCarloSim(float T, int L, int n, float p) : n{n}, time{0}, is{new ImpureIsingSubstance(T,L,p)}, distr{0, L-1}, pRatios{
            (Eigen::Array<float, 5, 1>{8,4,0,-4,-8}/(-T)).exp()
        } {};
        MonteCarloSim(std::shared_ptr<ImpureIsingSubstance> isptr, int n) : n{n}, time{0}, is{isptr}, distr{0, isptr->L-1}, pRatios{
            (Eigen::Array<float, 5, 1>{8,4,0,-4,-8}/(-isptr->T)).exp()
        } {};
        Eigen::ArrayXf magnetization{n};
        Eigen::ArrayXf energy{n};
        Eigen::ArrayXf autocorr{n};
        Eigen::ArrayXi timeArray = Eigen::ArrayXi::LinSpaced(n,0,n-1);
        float getHeatCapacity(){return (energy.segment(4*lag,time-4*lag)-energy.segment(4*lag,time-4*lag).mean()).matrix().squaredNorm()/(is->T*is->T);}
        int getImpurityIndices(int i){return is->getIndices(i);}
        void runSimulation();
    };

    class TemperatureRun{
        float startT;
        float diffT;
        int m;
        float endT = startT+m*diffT;
        int step;
        int L;
        int n;
        float p;
        bool impurities{false};
        void temperatureStep();
        void smoothStep();
        int maxIndex(Eigen::ArrayXf v);
        public:
        TemperatureRun(float startT, int m, float diffT = 0.02,  int L=40, int n=20000) : startT{startT}, diffT{diffT}, n{n}, L{L}, m{m}, step{0}, impurities{false} {}; //Constructor for no impurities
        TemperatureRun(float startT, int m, float p, float diffT = 0.02,  int L=40, int n=20000) : startT{startT}, diffT{diffT}, n{n}, L{L}, m{m}, step{0}, p{p}, impurities{true} {}; //Constructor for impurities
        void goForARun();
        void smoothRun(int k);
        Eigen::ArrayXf heatCapacities = Eigen::ArrayXf::Zero(m);
        Eigen::ArrayXf tempArray = Eigen::ArrayXf::LinSpaced(m,startT,endT);
        int getMaxHCIndex(){return maxIndex(heatCapacities);}
    };

    class SimulatedAnnealing{
        /* Transition probability for lattice is >exp(-8/T). For geometric annealing we want
        high probability of acceptance at our temperature, say 0.95 which gives 
        T =-8/ln(0.95)~156, 0.9 -> T = 75, or 0.8 -> T = 36. */
        float T; //Current temperature, starts "high"
        float Tstep; //Cooling dump in temperature, if linear cooling
        float alpha; // Between 0 and 1 (chosen between 0.8 and 0.99) geometric cooling factor
        int k; //Current number of iterations without state changing
        int stopCriterion; //Number of iterations without state changing needed for termination
        std::shared_ptr<ImpureIsingSubstance> is; // We want this substance to be the same i.e. keep state from
        // one simulation to the next.
        int n; // Number of sweeps per temperature
        void annealingStep();
        Eigen::ArrayXf energy = Eigen::ArrayXf::Zero(0);
        int count = 0;
        public:
        SimulatedAnnealing(float T0, int L, float alpha, int stopCriterion, float p, int n) : T{T0}, alpha{alpha}, k{0}, n{n}, stopCriterion{stopCriterion}, is{std::make_shared<ImpureIsingSubstance>(T0, L, p)} {};
        void runSimulatedAnnealing();
        std::vector<std::vector<double>> getLattice();
    };
}