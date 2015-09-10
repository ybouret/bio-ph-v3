-- -----------------------------------------------------------------------------
--
-- times
--
-- -----------------------------------------------------------------------------
ftol   = 1e-4; -- differential fractional tolerance
diff_h = 1e-5; -- initial adaptive time step between to time steps

dt      = 0.1;
dt_save = 0.2;
t_run   = 60+600+600;

T       = 298;

Em = 0;
zeta0 = Em * F / (R*T);

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
    {"HCO3-",-1},
    {"CO3--",-2},
    {"HY",    0},
    {"Y-",   -1},
    {"OSM",   0}
};

-- -----------------------------------------------------------------------------
-- Thermo
-- -----------------------------------------------------------------------------
Kw      = 1.0e-14;
K_Henry = 29.41;   -- atm/(mol/L)

K1      = 4.45e-7;
K2      = 5.6e-11;

pKY = 6.2;

function P_CO2(t)
local  P0    = 40.0;
--local  W     = 60;
--if (t>=5) and (t<=W*math.pi+5) then
--  return P0/760.0 + (40.0/760.0) * math.sin((t-5)/W)^2;
--end
return P0/760;
end

function kappa(t)
local  CO2_aq = P_CO2(t) / K_Henry;
return CO2_aq * K1;
end


-- -----------------------------------------------------------------------------
--
-- Equilibria
--
-- -----------------------------------------------------------------------------

eqs =
{
    { "water", 1e-14, { 1, "H+"}, {  1, "HO-"} },
    { "bicarb", "kappa",          {  1, "H+" }, { 1, "HCO3-"} },
    { "carb",   K2,               {  1, "H+" }, { 1, "CO3--"}, { -1, "HCO3-" } },
    { "buffer", 10^(-pKY),        { -1, "HY" }, { 1, "Y-"   }, {1, "H+" } },
};

eff =
{
    "lambda_Na",
    "lambda_K",
    "lambda_Cl",
    "NaK",
    "NHE",
    "AE2"
    -- , "inject"
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
    { "osmolarity", 300e-3},
    { 0, { 1, "HY"}, {1, "Y-" } }
};


-- -----------------------------------------------------------------------------
--
-- Outside
--
-- -----------------------------------------------------------------------------

init0 =
{
    { 10^(-7.4),    { 1, "H+"  } },
    { 4e-3,         { 1, "K+"  } },
    { 140e-3,       { 1, "Na+" } },
    { 100e-3,       { 1, "Cl-" } },
    { "osmolarity", 300e-3 },
    { 0, { 1, "HY"}, {1, "Y-" } }
}

init1 =
{
    { 10^(-6.5),    { 1, "H+"  } },
    { 4e-3,         { 1, "K+"  } },
    { 140e-3,       { 1, "Na+" } },
    { 100e-3,       { 1, "Cl-" } },
    { "osmolarity",  300e-3 },
    { 0, { 1, "HY"}, {1, "Y-" } }
}


out =
{
    "init0",
    "init1"
};


function weights(t)
--local  w = exp(-t);
--return w,1-w;
if (t<=60) then
return 1,0
end
if(t<=60+600)then
return 0,1
end
return 1,0
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

-- -----------------------------------------------------------------------------
-- Common surface capacitance
-- -----------------------------------------------------------------------------


Cm = 1; -- muFarad/cm^2

-- -----------------------------------------------------------------------------
-- Article from which the permeabilities were fitted
-- -----------------------------------------------------------------------------


Cm      = 1;   -- muF/cm^2

Capa_exp = 207e-12;               -- from article 207 pF
Surf_exp = Capa_exp / (1e-14*Cm); -- in micron^2

-- -----------------------------------------------------------------------------
--
-- Geometry
--
-- -----------------------------------------------------------------------------


a = 12.5; -- in microns
b = 5;    -- in microns
c = 5;    -- in microns

surface = EllipsoidSurface(a,b,c);
volume  = EllipsoidVolume(a,b,c);


-- -----------------------------------------------------------------------------
-- INWARD Na in  moles/s/m^2
-- -----------------------------------------------------------------------------
function lambda_Na(t,Cin,Cout,params)
local zeta = params["zeta"];
local zz   = zeta;
local Na   = "Na+";
local S    = params["S"];
local V    = params["V"];

a = {}
local Perm = SP_Na(t,zeta)/Surf_exp;                    -- in microns/s
local Flux = Perm * Psi(zz)*(Cout[Na]-Cin[Na]*exp(zz)); -- in moles/L*microns/s
local J    = 1e-3 * Flux;                               -- in moles/m^2/s
a[Na]      = J;
return a;

end

-- -----------------------------------------------------------------------------
-- INWARD chloride moles/s/m^2
-- -----------------------------------------------------------------------------
function lambda_Cl(t,Cin,Cout,params)
local zeta = params["zeta"];
local zz   = -zeta;
local Cl = "Cl-";

a = {};
local Perm = SP_Cl(t,zeta)/Surf_exp;                    -- in microns/s
local Flux = Perm * Psi(zz)*(Cout[Cl]-Cin[Cl]*exp(zz)); -- in moles/L*microns/s
local J    = 1e-3*Flux;                                 -- in moles/m^2/s
a[Cl]      = J;
return a;

end


-- -----------------------------------------------------------------------------
-- INWARD potassium moles/s/m^2
-- -----------------------------------------------------------------------------
function lambda_K(t,Cin,Cout,params)
local zeta = params["zeta"];
local zz   = zeta;
local K  = "K+";

a = {};
local Perm = SP_K(t,zeta)/Surf_exp;                    -- in microns/s
local Flux = Perm * Psi(zz)*(Cout[K]-Cin[K]*exp(zz));  -- in moles/L*microns/s
local J    = 1e-3*Flux;                                -- in moles/m^2/s
a[K] = J;
return a;

end

-- -----------------------------------------------------------------------------
-- NaK/ATPase INWARD moles/s, unscaled
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
-- NHE INWARD moles/s, unscaled
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


-- -----------------------------------------------------------------------------
-- AE2 INWARD moles/s
-- -----------------------------------------------------------------------------
K_AE = 10e-3;
function AE2(t,Cin,Cout,params)
local b   = Cin["HCO3-"];
local rho = b/(b+K_AE);
a = {}
a["Cl-"]   = rho;
a["HCO3-"] = -rho;
return a;
end

-- -----------------------------------------------------------------------------
-- Changing potential
-- -----------------------------------------------------------------------------

function inject(t,Cin,Cout,params)
a = {}
a["Na+"] = 0.01;
return a;
end




