 #include "IsingSubstance.h"
#include "MonteCarlo.h"

using namespace std;

namespace IsingProject {
    void MonteCarloSim::newState(){
        /*Generates random 2d-indices in lattice, calculates energy-difference by flipping spin/state change
        at chosen index according to 2d ising model, finds the probability ratio from pRatios array and uses 
        this in the Metropolis algorithm. For ImpureIsingSubstance, isItPure uses polymorphism to
        check whether the index (column-major) belongs to is->ImpureIndices. If so, fixed spin means
        no flip can occur.*/
        int* l_ptr = is->lattice.data();
        //const int L = is.L;
        for(int i=0; i<is->L*is->L; i++){
            int ix = distr(gen);
            int iy = distr(gen);
            int ind = ix + iy*is->L;
            if(!is->isItPure(ind)){continue;}
            //Finding change in energy:
            int count = l_ptr[((iy + 1)%is->L)*is->L + ix] + l_ptr[(ix - 1 + is->L)%is->L + is->L*iy] + l_ptr[ix+ is->L*((iy - 1 + is->L)%is->L)] + l_ptr[(ix + 1)%is->L + is->L*iy];
            //Probability ratio:
            int rIndex = static_cast<int>(2-0.5*l_ptr[ind]*count);
            //Metropolis algorithm:
            float pr = pRatios(rIndex);
            float draw = dis(gen);
            if(pr >= 1 || pr > draw){
                //Flipping spin according to Metropolis:
                l_ptr[ind]*=-1;
            }

        }
    }

    void MonteCarloSim::timeStep(){
        newState();
        magnetization(time)=is->magnetizationCalc();
        energy(time)=is->energyCalc();
        autocorr(time)=std::abs(autoCorrolation());
        time++;
    }

    float MonteCarloSim::autoCorrolation(){
        if(time > 3*lag){
            float mean = magnetization.segment(lag, time-lag).mean();
            auto mag = (magnetization.segment(2*lag,time-2*lag)-mean).matrix();
            auto mag_delay = (magnetization.segment(lag, time-2*lag)-mean).matrix();
            float cov = mag.dot(mag_delay);
            float var = ((magnetization.segment(lag, time-lag) - mean).matrix()).squaredNorm();
            return cov/var*(time-lag)/(time-2*lag);
        }
        else{return 0;}
    }

    void MonteCarloSim::runSimulation(){
        time=0;
        while(time < n){
            timeStep();
        }
    }

    void TemperatureRun::temperatureStep(){
        MonteCarloSim sim = impurities ? MonteCarloSim(startT+diffT*step,L,n, p) : MonteCarloSim(startT+diffT*step,L,n);
        sim.runSimulation();
        heatCapacities(step)= sim.getHeatCapacity();
        step++;
    }

    void TemperatureRun::goForARun(){
        while(step<m){temperatureStep(); cout << step << endl;}
    }

    void TemperatureRun::smoothStep(){
        MonteCarloSim sim = impurities ? MonteCarloSim(startT+diffT*step,L,n, p) : MonteCarloSim(startT+diffT*step,L,n);
        sim.runSimulation();
        heatCapacities(step) +=sim.getHeatCapacity();
        std::cout << heatCapacities(step) << std::endl;
        step++;
    }

    void TemperatureRun::smoothRun(int k){
        for(int i = 0; i < k; i++){
            while(step<m){smoothStep();}
            step = 0;
        }
        heatCapacities/=k;
    }

    int TemperatureRun::maxIndex(Eigen::ArrayXf v){
        int ind = 0;
        int m = v[0];
        for(int i = 1; i < v.size(); i++){if(v[i] > m){m=v[i];ind = i;}}
        return ind;
    }

    void SimulatedAnnealing::annealingStep(){
        Eigen::ArrayXXi initLattice = is -> lattice; //Testing whether the lattice changes each step
        is->T=T;
        MonteCarloSim mcs{is, n};
        //cout << " " << mcs.getImpurityIndices(500);
        mcs.runSimulation();
        T*=alpha;
        if((initLattice == is -> lattice).all()){
            k+= 1;
            cout << " hello ";
        }
        else{k=0;}
        auto diff_mask = initLattice != is->lattice;
        for (int i = 0; i < diff_mask.size(); ++i) {
            if (diff_mask(i)) {cout << i << endl;}
            }
        cout << mcs.pRatios(0) << " " << mcs.pRatios(1) << " " << mcs.pRatios(2) << " " << mcs.pRatios(3) << " " << mcs.pRatios(4) << endl;
    }

    void SimulatedAnnealing::runSimulatedAnnealing(){
        int c = 0;
        while(k<stopCriterion && c < 150){
            annealingStep();
            cout << "__" << c << "__" << endl;
            c++; //they said the thing!
        }
        //for(auto obj : is->impurityIndices){
        //    cout << obj << endl;
        //}
    }

    std::vector<std::vector<double>> SimulatedAnnealing::getLattice(){
        std::vector<std::vector<double>> m(is->L, std::vector<double>(is->L,0));
        for(int i = 0; i<is->L; i++){
            for(int j = 0; j<is->L; j++){
                m.at(i).at(j) = static_cast<double>(is->lattice(i,j));
            }
        }
        return m;
    }
}
