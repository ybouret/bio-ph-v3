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


void Cell:: compute_fields( array<double> &dYdt, double t, const array<double> &Y )
{
    load_state(t, Y);
    const double zeta = (Em*__Faraday__)/(__R__*Temperature);
    
    // chemical part
    compute_rates(t,zeta);
    sol_tmp->save(dYdt);
    if(false)
    {
        eqs.load_C(Y);
        eqs.load_dC(dYdt);
        
        eqs.legalize_dC(t);
        eqs.save_dC(dYdt);
    }
    // physical part
    dYdt[idxE] = 0;
    
}

void Cell:: step( double t1, double t2 )
{
    //-- prepare initial state
    sol_ins->save(X);
    X[idxE] = Em;
    compute_out(t1);
    
    //-- internal computation
    odeint(drvs,X,t1,t2,ctrl,NULL);
    
    
    //-- save final state
    sol_ins->load(X);
    Em = X[idxE];
    compute_out(t2);
    
}
