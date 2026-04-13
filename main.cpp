#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "IsingSubstance.h"
#include "MonteCarlo.h"
#include <chrono>
#include <matplot/matplot.h>
using namespace std;
using namespace IsingProject;

int main()
{
   auto start = std::chrono::steady_clock::now();
   /*
   //Code for generating timeseries of energy and magnetization
   constexpr int L = 40;
   constexpr float T = 1.7f;
   IsingProject::MonteCarloSim isingSim(T,L,200000);
   isingSim.runSimulation();

   std::cout << "Heat Capacity:" << isingSim.getHeatCapacity() << std::endl;

   matplot::plot(isingSim.timeArray, isingSim.autocorr);
   matplot::show();

   matplot::plot(isingSim.timeArray, isingSim.energy);
   matplot::hold(true);
   matplot::plot(isingSim.timeArray, isingSim.magnetization);

   matplot::ylim({-2.0,1.0});
   matplot::show();
   */

   /* //Code for plotting heat capacity and calculating max C.
   int L = 40;
   float startTemp = 2.1f;
   float tempStep = 0.02f;
   float endTemp = 2.5f;
   int m = static_cast<int>((endTemp-startTemp)/tempStep);
   int n = 5000;

   L=10;
   TemperatureRun tr10{startTemp, m, 0.02, L, n};
   tr10.smoothRun(10);
   L=20;
   TemperatureRun tr20{startTemp, m, 0.02, L, n};
   tr20.smoothRun(20);
   L=40;
   TemperatureRun tr40{startTemp, m, 0.02, L, n};
   tr40.smoothRun(30);
   L=100;
   TemperatureRun tr100{startTemp, m, 0.02, L, n};
   tr100.smoothRun(30);
   cout << tr10.getMaxHCIndex() << " " <<tr10.getMaxHCIndex()*tempStep+startTemp << " " << tr10.heatCapacities[tr10.getMaxHCIndex()] << endl;
   cout << tr20.getMaxHCIndex() << " " <<tr20.getMaxHCIndex()*tempStep+startTemp << " " << tr20.heatCapacities[tr20.getMaxHCIndex()] << endl;
   cout << tr40.getMaxHCIndex() << " " <<tr40.getMaxHCIndex()*tempStep+startTemp << " " << tr40.heatCapacities[tr40.getMaxHCIndex()] << endl;
   cout << tr100.getMaxHCIndex() << " " <<tr100.getMaxHCIndex()*tempStep+startTemp << " " << tr100.heatCapacities[tr100.getMaxHCIndex()] << endl;

   std::cout << "Time: " << elapsed.count() << "s\n";

   matplot::plot(tr10.tempArray, tr10.heatCapacities);
   matplot::hold(true);

   matplot::plot(tr20.tempArray, tr20.heatCapacities);
   matplot::hold(true);
   
   matplot::plot(tr40.tempArray, tr40.heatCapacities);
   matplot::hold(true);

   matplot::plot(tr100.tempArray, tr100.heatCapacities);
   matplot::show();
   matplot::hold(false);

   matplot::semilogy(tr10.tempArray, tr10.heatCapacities);
   matplot::hold(true);

   matplot::semilogy(tr20.tempArray, tr20.heatCapacities);
   matplot::hold(true);
   
   matplot::semilogy(tr40.tempArray, tr40.heatCapacities);
   matplot::hold(true);
   matplot::semilogy(tr100.tempArray, tr100.heatCapacities);
   matplot::show();
   matplot::hold(false);

   matplot::plot(tr10.tempArray, tr10.heatCapacities.array().log().eval());
   matplot::hold(true);
   matplot::plot(tr20.tempArray, tr20.heatCapacities.array().log().eval());
   matplot::hold(true);
   matplot::plot(tr40.tempArray, tr40.heatCapacities.array().log().eval());
   matplot::hold(true);
   matplot::plot(tr100.tempArray, tr100.heatCapacities.array().log().eval());
   matplot::show();
   */
   /*
   constexpr int L = 15;
   constexpr float T = 2.3f;
   int n = 10000;
   float p = 0.02;

   IsingProject::MonteCarloSim isingSim(T,L,n,p);
   isingSim.runSimulation();

   std::cout << "Heat Capacity:" << isingSim.getHeatCapacity() << std::endl;
   isingSim.printImpurityIndices();

   matplot::plot(isingSim.timeArray, isingSim.autocorr);
   matplot::show();

   matplot::plot(isingSim.timeArray, isingSim.energy);
   matplot::hold(true);
   matplot::plot(isingSim.timeArray, isingSim.magnetization);

   matplot::ylim({-2.0,1.0});
   matplot::show();
   matplot::hold(false);
   */
   /* //Simulated ising lattice with impurities percentage p
   float startTemp = 2.1f;
   float tempStep = 0.01f;
   float endTemp = 2.5f;
   int m = static_cast<int>((endTemp-startTemp)/tempStep);
   n = 5000;

   TemperatureRun trImp(startTemp, m, p, tempStep, L, n);
   trImp.smoothRun(70);
   matplot::hold(true);
   matplot::plot(trImp.tempArray, trImp.heatCapacities);

   TemperatureRun tr(startTemp, m, tempStep, L, n);
   tr.smoothRun(70);

   matplot::plot(tr.tempArray, tr.heatCapacities);

   
   matplot::ylim({0,tr.heatCapacities.maxCoeff() + 1});
   matplot::show();


   */

   float T0 = 750; int L = 40; float alpha = 0.9; float p = 0.40; int n = 5000;
   SimulatedAnnealing sa{T0, L, alpha, 5, p, n};
   sa.runSimulatedAnnealing();
   
   matplot::image(sa.getLattice());
   matplot::show(); 
   
   auto end = std::chrono::steady_clock::now();
   std::chrono::duration<double> elapsed = end - start;

   return 0;
} 
