#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "luavm.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/memory/pooled.hpp"

class HVariables : public object
{
public:
    static const char  *NAMES_REG[]; //!< built-in parameter names: zeta, V, S
    static const char  *LOADS_REG[]; //!< built-in parameter loads: zeta0, volume, surface
    static const size_t EXTRA_NUM;   //!< the number

    typedef vector<string,memory::pooled::allocator> vector_s;

    vector_s names;
    vector_s loads;
    
    explicit HVariables();
    virtual ~HVariables() throw();

protected:
    void ReleaseContent() throw();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HVariables);
};

//! prototype minimal cell
/**
 - concentrations:      mol/L
 - volume:              mu^3 = 1e-18 m^2 = 1e-15 L
 - surface:             mu^2
 - potential:           mV (but use zeta instead)
 - surface capacitance: muF/cm^2, 1e-14 F/mu^2
 - Fluxes: mol/m^2/s
 */
class HCell : public HVariables
{
public:
    static const double ZETA_MAX;           //!< max value of acceptable zeta
    
    //!
    /**
     - Loading species from           "lib"
     - Loading equilibria from        "eqs"
     - Loading effectors from         "eff"
     - Loading inside0   from         "ini"
     - Loading outside solutions from "out"
     - The global initial volume and surface must be available
     */
    explicit HCell(lua_State   *vm,
                   const double t0,
                   const char  *extra_params_reg[],
                   const size_t extra_params_num);

    virtual ~HCell() throw();

    lua_State         *L;           //!< internal virtual machine

    //__________________________________________________________________________
    //
    // Chemical Part
    //__________________________________________________________________________
    __lua::Library     lib;         //!< the library
    __lua::Equilibria  eqs;         //!< the global chemical system
    const size_t      &N;           //!< #eqs
    const size_t      &M;           //!< #species
    variables          params;      //!< extra parameters, initially 0
    const size_t      &nvar;        //!< params.size
    const size_t       iZeta;       //!<
    const size_t       iVolume;     //!<
    const size_t       iSurface;    //!< 
    __lua::Effectors   eff;         //!< effectors, must FLUXES !!!
    vector_t           inside;      //!< initial inside concentration (+extra vars)
    double             tmx;         //!< current time
    vector_t           rho;         //!< temporary rates, #M
    vector_t           in;          //!< current inside
    matrix_t           outside;     //!< possible outside solutions
    vector_t           out;         //!< resulting from mix, using the lua "weights" function
    vector_t           weights;     //!< to store weights of outsides compositions

    //__________________________________________________________________________
    //
    // Physical Part
    //__________________________________________________________________________
    const double       T;           //!< temperature
    const double       E2Z;         //!< zeta  = E2Z * Em, E2Z=F/(R*T)
    const double       Z2E;         //!< Em    = Z2E * zeta
    const double       Cm;          //!< surface capacitance

    //__________________________________________________________________________
    //
    // differential part
    //__________________________________________________________________________
    diff_solver   odeint; //!< solver, initialized for nvar
    const double  diff_h; //!< initial time step for each dt
    diff_equation diffeq; //!< use Call, calling virtual Rates...
    size_t        ncalls; //!< internal counter

    //__________________________________________________________________________
    //
    // control functions
    //__________________________________________________________________________

    //! using outside solutions...
    /**
     The result is set in out
     */
    void   ComputeOutsideComposition(const double t);



    //! compute no-flux zeta
    double ComputeRestingZeta(const double t);


    //! compute all chemical rates
    /**
     the change in zeta is automatically computed
     */
    virtual void Rates( array<double> &dYdt, double t, const array<double> &Y ) = 0;

    //! solve the Rates equation from t0 to t1
    void Step(array<double> &Y, double t0, double t1);

    //__________________________________________________________________________
    //
    // I/O functions
    //__________________________________________________________________________
    ios::ostream & add_header( ios::ostream &fp ) const;
    ios::ostream & add_values( ios::ostream &fp, const array<double> &Y ) const;

    ios::ostream & add_header_csv( ios::ostream &fp ) const;
    ios::ostream & add_values_csv( ios::ostream &fp, const array<double> &Y ) const;


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
    const array<string> & fill_params_reg(const char  *extra_params_reg[],
                                          const size_t extra_params_num);

    void Call( array<double> &dYdt, double t, const array<double> &Y );

    //! compute volumic charge rate
    /**
     assuming out is computed, and using in
     \return the algebraic signed flux, moles/s/m^2
     */
    double ComputeVolumicChargeRate(double zeta);

};


#endif

