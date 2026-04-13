#include "IsingSubstance.h"

namespace IsingProject {
    Eigen::ArrayXXi IsingSubstance::setSumOfNeighbors(){
        Eigen::ArrayXXi shiftLeft{L,L};
        Eigen::ArrayXXi shiftUp{L,L};
        Eigen::ArrayXXi shiftRight{L,L};
        Eigen::ArrayXXi shiftDown{L,L};

        shiftUp.topRows(L-1)=lattice.block(L-1,L,1,0);
        shiftLeft.leftCols(L-1)=lattice.block(L,L-1,0,1);
        shiftRight.rightCols(L-1)=lattice.block(L,L-1,0,0);
        shiftDown.bottomRows(L-1)=lattice.block(L-1,L,0,0);

        shiftUp.bottomRows(1)=lattice.row(0);
        shiftLeft.rightCols(1)=lattice.col(0);
        shiftRight.leftCols(1)=lattice.col(L);
        shiftDown.topRows(1)=lattice.row(L);
        
        return shiftDown+shiftUp+shiftLeft+shiftRight;
    };

    float IsingSubstance::energyCalc(){
        /*Outputs energy as a scalar of current grid.*/
        Eigen::ArrayXXi shiftLeft{L,L};
        Eigen::ArrayXXi shiftUp{L,L};
        shiftUp.topRows(L-1)=lattice.block(1,0,L-1,L);
        shiftLeft.leftCols(L-1)=lattice.block(0,1,L,L-1);

        shiftUp.bottomRows(1)=lattice.row(0);
        shiftLeft.rightCols(1)=lattice.col(0);


        return -(lattice*(shiftLeft+shiftUp)).sum()*Ninv;
    };

    float IsingSubstance::magnetizationCalc(){
        return std::abs(lattice.sum()*Ninv);
    }

    void ImpureIsingSubstance::setImpurities(){
        int i = 0;
        std::vector<int> taken;
        std::random_device rd; 
        std::mt19937 gen{rd()}; 
        std::uniform_int_distribution<> distr{0, L*L-1}; 
        while(i < impurities){
            int newindex = distr(gen);
            if(std::find(taken.begin(), taken.end(), newindex) == taken.end()){
                taken.push_back(newindex);
                impurityIndices[i] = newindex;
                i++;
            }
        }
    }
}
