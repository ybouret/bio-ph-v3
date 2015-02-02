#ifndef H_CELL_INCLUDED
#define H_CELL_INCLUDED 1

#include "luavm.hpp"
#include "yocto/physics/constants.hpp"
#include "yocto/ios/ostream.hpp"

class HCell
{
public:
    static const char  *PARAMETERS[];
    static const size_t NUM_PARAMS;
    static const double ZETA_MAX;


    //!
    /**
     - Loading species from           "lib"
     - Loading equilibria from        "eqs"
     - Loading effectors from         "eff"
     - Loading inside0   from         "ini"
     - Loading outside solutions from "out"
     - The global initial volume and surface must be available
     */
    explicit HCell(lua_State    *vm,
                   const double  t0);
    virtual ~HCell() throw();

    lua_State        *L;
    __lua::Library    lib;         //!< the library
    __lua::Equilibria eqs;         //!< the global chemical system
    const size_t      N;           //!< #eqs
    const size_t      M;           //!< #species
    parameters        params;      //!< extra parameters
    const size_t      iZeta;       //!< zeta index
    const size_t      iVolume;     //!< volume index
    const size_t      iSurface;    //!< surface index
    const size_t      iActiveS;    //!< active surface, initially surface
    const size_t      nvar;        //!< #nvar for all
    __lua::Effectors  eff;         //!< effectors
    vector_t          inside;      //!< initial inside concentration (+extra vars)
    matrix_t          outside;     //!< possible outside solutions
    vector_t          out;         //!< resulting from mix, using the lua "weights" function
    vector_t          weights;     //!< to store wetighs
    vector_t          in;          //!< temporary, nvar, initialized to inside0
    double            tmx;         //!< local time
    vector_t          rho;         //!< temporary rates, #M
    const double      Temperature; //!<
    const double      E2Z;         //!< F/(RT)
    const double      Z2E;         //!< (RT)/F;

    diff_solver   odeint; //!< solver, initialized for nvar
    const double  diff_h; //!< initial time step for each dt
    diff_equation diffeq; //!< use Call, calling virtual Rates...
    size_t        ncalls; //!< internal counter

    //! using outside solutions...
    /**
     The result is set in out
     */
    void   ComputeOutsideComposition(const double t);

    //! using in and out
    /**
     store the result in inside[iZeta] as well
     */
    double ComputeRestingZeta(const double t);

    //! to be implemented
    virtual void   Rates( array<double> &dYdt, double t, const array<double> &Y ) = 0;

    //! solve the Rates equation from t0 to t1
    void Step(array<double> &Y, double t0, double t1);

    void add_header( ios::ostream &fp ) const;
    void add_values( ios::ostream &fp, const array<double> &Y ) const;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(HCell);
    //! compute individual fluxes from effectors
    /**
     \return the algebraic signed flux
     */
    double ComputeFluxes(double zeta);
    void   Call( array<double> &dYdt, double t, const array<double> &Y );

};


#endif
