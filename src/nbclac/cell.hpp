#ifndef PLOS_CELL_INCLUDED
#define PLOS_CELL_INCLUDED 1

#include "../h+core/hcell.hpp"
#include "yocto/lua/lua-maths.hpp"

class Cell : public HCell
{
public:
    const size_t iNa;
    const size_t iK;
    const size_t iCl;
    const size_t iH;
    const size_t iB;
    const size_t iLacH;
    
    effector &leak_K;
    effector &leak_Na;
    effector &leak_Cl;

    effector     &NHE;
    effector     &AE2;
    effector     &NaK;
    const double  alpha;    //!< NHE fraction in [0:1]
    effector     &NBC;
    effector     &MCT1;
    effector     &MCT4;
    const double  Vm1; //!< max speed for MCT1
    const double  Vm4; //!< max speed for MCT4

    Lua::Function<double> Lambda;

    explicit Cell( lua_State *vm, const double t);
    virtual ~Cell() throw();

    //! evaluate SteadyState Zeta @t=0
    /**
     use channels/NaK using in=inside0 and computed out @t=0
     */
    double SteadyStateZeta();

    //! set Potential to Em using lambda_K.pace @t=0
    void   SetSteadyStatePotential(double Em);

    //! SetSteady state and all paces
    void Setup(double Em);

    void Powers(double t, const array<double> &Y, double &p1, double &p4);


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);

    //! 1.5*lambda_K + lambda_Na - lambda_Cl-
    double ComputePassiveZeroFlux(double zeta);
    double ComputePassiveZeroKFlux(double pace);

    virtual void   Rates( array<double> &dYdt, double t, const array<double> &Y );

    
};

#endif

