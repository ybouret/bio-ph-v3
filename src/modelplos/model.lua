-- -----------------------------------------------------------------------------
--
-- standard
--
-- -----------------------------------------------------------------------------
exp  = math.exp;
tanh = math.tanh;


-- -----------------------------------------------------------------------------
--
-- Components
--
-- -----------------------------------------------------------------------------

lib =
{
    {"H+",   1},
    {"HO-", -1},
    {"Na+",  1},
    {"Cl-", -1},
    {"K+",   1},
    {"OSM",  0}
};


-- -----------------------------------------------------------------------------
--
-- Equilibria
--
-- -----------------------------------------------------------------------------

eqs =
{
    { "water", 1e-14, { 1, "H+"}, { 1, "HO-"} }
};

eff =
{
    "lambda_Na",
    "lambda_K",
    "lambda_Cl",
    "NaK"
};


-- -----------------------------------------------------------------------------
--
-- Inside
--
-- -----------------------------------------------------------------------------

ini =
{
    { 10^(-7.2), { 1, "H+" } },
    { 140e-3,    { 1, "K+"  } },
    { 10e-3,     { 1, "Na+" } },
    { 20e-3,     { 1, "Cl-" } },
    { "osmolarity", 300e-3}
};


-- -----------------------------------------------------------------------------
--
-- Outside
--
-- -----------------------------------------------------------------------------

init0 =
{
    { 10^(-7.4),    { 1, "H+"  } },
    -- { 0,            { 1, "HY"  }, { 1, "Y-" } },
    { 4e-3,         { 1, "K+"  } },
    { 140e-3,       { 1, "Na+" } },
    { 100e-3,       { 1, "Cl-" } },
    { "osmolarity", 300e-3 }
}

init1 =
{
    { 10^(-6.5),    { 1, "H+"  } },
    --{ 0,            { 1, "HY"  }, { 1, "Y-" } },
    { 4e-3,         { 1, "K+"  } },
    { 140e-3,       { 1, "Na+" } },
    { 100e-3,       { 1, "Cl-" } },
    { "osmolarity",  500e-3 }
}


out =
{
    "init0",
    "init1"
};


function weights(t)
return 0.8,0.2
end

-- -----------------------------------------------------------------------------
-- Helper
-- -----------------------------------------------------------------------------

function TableSumValues(X)
ans = 0;
for key,value in pairs(X) do ans = ans+value; end
return ans;
end

-- -----------------------------------------------------------------------------
-- Surface * Permeability in micron^3 / s
-- -----------------------------------------------------------------------------

function SP_K( t, x )
return ((3.70714279) * exp( (0.4310085034) * x ) + (1.507274188));
end

-- -----------------------------------------------------------------------------
-- Surface * Permeability in micron^3 / s
-- -----------------------------------------------------------------------------
function SP_Cl(t,x)
return (1.194531961) * exp( (0.05658342092) * x );
end

-- -----------------------------------------------------------------------------
-- Surface * Permeability in micron^3 / s
-- -----------------------------------------------------------------------------

function SP_Na(t,x)
return (1.360058333) * exp( (0.05223441992) * x );
end

Cm = 10 * 1.0e-15; -- Farad/ micron^2

Capa_exp = 207e-12;       -- from article 207 pF
Surf_exp = Capa_exp / Cm; -- in micron^2

-- -----------------------------------------------------------------------------
--
-- Geometry
--
-- -----------------------------------------------------------------------------


a = 12.5;
b = 5;
c = 5;

surface = EllipsoidSurface(a,b,c);
volume  = EllipsoidVolume(a,b,c);

-- channels extrapolations
passive_ratio = surface/Surf_exp;

-- -----------------------------------------------------------------------------
-- INWARD Na moles/s
-- -----------------------------------------------------------------------------
function lambda_Na(t,Cin,Cout,params)
local zeta = params["zeta"];
local zz   = zeta;
local Na = "Na+";
a = {};
local rho = Psi(zz)*(Cout[Na]-Cin[Na]*exp(zz)) * SP_Na(t,zeta);
a[Na] = rho*passive_ratio;
return a;
end

-- -----------------------------------------------------------------------------
-- INWARD potassium moles/s
-- -----------------------------------------------------------------------------
function lambda_K(t,Cin,Cout,params)
local zeta = params["zeta"];
local zz   = zeta;
local K  = "K+";
a = {};
local rho = Psi(zz)*(Cout[K]-Cin[K]*exp(zz)) * SP_K(t,zeta);
a[K] = rho*passive_ratio;
return a;
end

-- -----------------------------------------------------------------------------
-- INWARD chloride moles/s
-- -----------------------------------------------------------------------------
function lambda_Cl(t,Cin,Cout,params)
local zeta = params["zeta"];
local zz   = -zeta;
local Cl = "Cl-";
a = {};
local rho = Psi(zz)*(Cout[Cl]-Cin[Cl]*exp(zz)) * SP_Cl(t,zeta);
a[Cl] = rho*passive_ratio;
return a;
end

-- -----------------------------------------------------------------------------
-- NaK/ATPase INWARD moles/s
-- -----------------------------------------------------------------------------
K_NaK = 12e-3;

function NaK(t,Cin,Cout,params)
local zeta  = params["zeta"];
local CNa   = Cin["Na+"];
local sig   = (CNa/(K_NaK+CNa)) * (1+tanh(0.39*zeta+1.28))*0.5;
local rho   = sig;
a = {}
a["K+"]  =  2*rho;
a["Na+"] = -3*rho;
return a;
end

-- -----------------------------------------------------------------------------
-- NHE INWARD moles/s
-- -----------------------------------------------------------------------------

L0    = 1000;
Kr    = 1.8e-8;
Kt    = 3.6e-6;
KNae  = 31e-3;
KHout = 2e-7;
function NHE(t,Cin,Cout,params)
local C    = Kr/Kt;
local x    = Cin["H+"]/Kr;
local xp1  = 1+x;
local cxp1 = 1+C*x;
local sig_num = x*xp1 + L0*C*x*cxp1;
local sig_den = L0*cxp1^2 + xp1^1;
local sig     = sig_num/sig_den;
local Nae     = Cout["Na+"];
local KNaeEff = KNae * (1+Cout["H+"]/KHout);
local sig_out = Nae / (KNaeEff+Nae);
sig = sig * sig_out;
ans = {}
ans["H+"]  = -sig;
ans["Na+"] =  sig;
return ans;
end




