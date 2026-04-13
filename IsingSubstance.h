#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <algorithm>

namespace IsingProject {
    class IsingSubstance{
        float Ninv;
        Eigen::ArrayXXi setSumOfNeighbors();
        public:
        float T;
        int L;
        float energyCalc();
        float magnetizationCalc();
        Eigen::ArrayXXi lattice;
        Eigen::ArrayXXi sumOfNeighbors{L,L};
        virtual bool isItPure(int i){return true;}
        virtual int getIndices(int i){return -1;}
        virtual int getIndSize(){return -1;}
        virtual void printImpurities(){};
        /*Possible eneIsingSubstance(float T, int L){energy-changes over differences in state, negative
        exponentiated with boltzmann factor (unit constants)
        to get ratio of probablities in Metropolis algoritm.*/
        IsingSubstance(float T, int L) : T{T}, L{L}, Ninv{ 1/static_cast<float>(L*L)}, lattice{
            (2*(Eigen::ArrayXXf::Random(L,L).floor()+0.5)).cast<int>()
        }
        {};
        virtual ~IsingSubstance() {};

    };

    class ImpureIsingSubstance : public IsingSubstance {
        float p;
        int impurities;
        void setImpurities();
        public:
        bool isItPure (int i) override {return !(i == impurityIndices).any();}
        int getIndices(int i){return impurityIndices(i);}
        void printImpurities(){std::cout << impurities;}
        int getIndSize(){return impurityIndices.size();}
        ImpureIsingSubstance(float T, int L, float prob) : IsingSubstance(T, L), p{prob}, impurities{static_cast<int>(std::floor(p*L*L))} {setImpurities();};
        Eigen::ArrayXi impurityIndices{impurities}; //1D column major lattice indices containing indices of the impurities
    };
}