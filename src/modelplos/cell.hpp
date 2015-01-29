#ifndef PLOS_CELL_INCLUDED
#define PLOS_CELL_INCLUDED 1

#include "../h+cell/cell.hpp"

class Cell : public HCell
{
public:
    const size_t iNa;
    const size_t iK;
    const size_t iCl;
    const size_t iH;

    explicit Cell( lua_State *vm, const double t);
    virtual ~Cell() throw();

    //! evaluate SteadyState Zeta @t=0
    /**
     use channels/NaK using in=inside0 and computed out @t=0
     */
    double SteadyStateZeta();

    //! set Potential to Em using lambda_K.pace @t=0
    void   SetSteadyStatePotential(double Em);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);

    //! 1.5*lambda_K + lambda_Na - lambda_Cl-
    double ComputePassiveZeroFlux(double zeta);
    double ComputePassiveZeroKFlux(double pace);

};

#endif

