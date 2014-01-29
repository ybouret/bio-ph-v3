#include "cell.hpp"


void Cell:: compute_rates(double t, double zeta)
{
    chemical::solution       &rates = *sol_tmp;
    const chemical::solution &S     = *sol_ins;
    const chemical::solution &S_out = *sol_out;
    
    leak(rates, t, zeta, S, S_out);
    append(rates, t, zeta, S, S_out);
}


void Cell:: save_state( array<double> &Y ) const
{
    assert( Y.size() == nvar );
    sol_ins->save(Y);
    Y[idxE] = Em;
}

void Cell:: load_state( double t, const array<double> &Y )
{
    assert( Y.size() == nvar );
    sol_ins->load(Y);
    Em = Y[idxE];
    
    compute_out(t);
    
}

void Cell:: NormalizeState( array<double> &Y, double t)
{
    eqs.load_C(Y);
    if(!eqs.normalize_C(t))
        throw exception("Invalid Concentrations");
    eqs.save_C(Y);
}


void Cell:: ComputeFields( array<double> &dYdt, double t, const array<double> &Y )
{
    load_state(t, Y);
    const double zeta = (Em*__Faraday__)/(__R__*Temperature);
    
    // collect chemical part
    compute_rates(t,zeta);
    sol_tmp->save(dYdt);
    
    // legalize
    {
        eqs.load_C(Y);
        eqs.load_dC(dYdt);
        
        eqs.legalize_dC(t);
        eqs.save_dC(dYdt);
        sol_tmp->load(dYdt);
    }

    
    // physical part
    const double zC    = sol_tmp->sum_zC();           // mol/L
    const double dCdt  = 1e-15 * zC;                  // mol/microns^3
    const double dQdt  = __Faraday__ * volume * dCdt; // net charges diff
    dYdt[idxE] = dQdt / Capa;
    
}

void Cell:: step( double t1, double t2 )
{
    //-- prepare initial state
    sol_ins->save(X);
    X[idxE] = Em;
    compute_out(t1);
    
    //-- internal computation
    odeint(drvs,X,t1,t2,ctrl,&cb);
    
    
    //-- save final state
    sol_ins->load(X);
    Em = X[idxE];
    compute_out(t2);
    
}
