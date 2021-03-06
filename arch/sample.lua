

local exp = math.exp;


-- -----------------------------------------------------------------------------
-- Species
-- -----------------------------------------------------------------------------
species =
{
    { "H+",     1, "SP_H" },
    { "HO-",   -1 },
    { "K+",     1 , "SP_K"  },
    { "Na+",    1,  "SP_Na" },
    { "Cl-",   -1,  "SP_Cl" },
    { "HCO3-", -1 },
    { "CO3--", -2 },
    { "HY",     0 },
    { "Y-",    -1 }
};

-- -----------------------------------------------------------------------------
-- Thermo
-- -----------------------------------------------------------------------------
Kw      = 1.0e-14;
K_Henry = 29.41;   -- atm/(mol/L)

K1      = 4.45e-7;
K2      = 5.6e-11;


function P_CO2(t)
local  P0    = 40.0;
--local  W     = 10;
--if (t>=0) and (t<=W*math.pi) then
--  return P0/760.0 + (40.0/760.0) * math.sin(t/W)^2;
--end
  return P0/760;
end

function kappa(t)
local  CO2_aq = P_CO2(t) / K_Henry;
return CO2_aq * K1;
end

pKY    = 6.2



equilibria =
{
    { "water",  Kw ,               {  1, "H+" }, { 1, "HO-"}   },
    { "bicarb", "kappa",           {  1, "H+" }, { 1, "HCO3-"} },
    { "carb",   K2,                {  1, "H+" }, { 1, "CO3--"}, { -1, "HCO3-" } },
    { "buffer",  10^(-pKY),        { -1, "HY" }, { 1, "Y-" }, {1, "H+" } },
}

-- ---------------------------------------------------------------------
-- Reduced Permeabilities, argument is zeta=FV/RT
-- ---------------------------------------------------------------------

-- Surface * Permeability in micron^3 / s
function SP_K( t, x )
return ((3.70714279) * exp( (0.4310085034) * x ) + (1.507274188));
end

-- Surface * Permeability in micron^3 / s
function SP_Cl(t,x)
return (1.194531961) * exp( (0.05658342092) * x );
end

-- Surface * Permeability in micron^3 / s

function SP_Na(t,x)
return (1.360058333) * exp( (0.05223441992) * x );
end

function SP_H(t,zeta)
return 0;
end

-- -----------------------------------------------------------------------------
-- Electric Parameters
-- -----------------------------------------------------------------------------
Cm = 10 * 1.0e-15; -- Farad/ micron^2

-- -----------------------------------------------------------------------------
-- Geometric Parameters
-- -----------------------------------------------------------------------------
surface,volume = SurfaceAndVolume(12.5,5,5);



-- -----------------------------------------------------------------------------
-- Internal Solution
-- -----------------------------------------------------------------------------


inside =
{
    { 10^(-7.2) ,  { 1, "H+" } },
    { 140e-3,      { 1, "K+"  } },
	{ 10e-3,       { 1, "Na+" } },
	{ 20e-3,       { 1, "Cl-" } },
    { 60e-3,       { 1, "HY"}, {1,"Y-"} }
};

-- ---------------------------------------------------------------------
-- External Solution(s)
-- ---------------------------------------------------------------------

init0 =
{
	{ 10^(-7.4),    { 1, "H+"  } },
	{ 0,            { 1, "HY"  }, { 1, "Y-" } },
	{ 4e-3,         { 1, "K+"  } },
	{ 140e-3,       { 1, "Na+" } },
	{ 100e-3,       { 1, "Cl-" } }
}

init1 =
{
	{ 10^(-6.5),    { 1, "H+"  } },
	{ 0,            { 1, "HY"  }, { 1, "Y-" } },
	{ 4e-3,         { 1, "K+"  } },
	{ 140e-3,       { 1, "Na+" } },
	{ 100e-3,       { 1, "Cl-" } }
}


outside =
{
    "init0",
    "init1"
};

function weights(t)
local t1 = 10
local t2 = 30*60+t1;
local t3 = 10+t2;

if(t<=0) then
return 1,0;
end

if(t<=t1) then
local w=t/t1;
return 1-w,w
end

if(t<=t2)then
return 0,1;
end

if(t<=t3)then
local w = (t-t2)/(t3-t2);
return w,1-w;
end

return 1,0;

end

-- -----------------------------------------------------------------------------
-- Effectors
-- -----------------------------------------------------------------------------
effectors = {
    "NaK",
    "AE",
    "NHE"
};


-- one function for effector
K_NaK = 12e-3;
function NaK(t,zeta,S,S_out)
local Na = S["Na+"];
local rho = 0.5*( 1+ math.tanh(0.39*zeta+1.28)) * (Na/(K_NaK+Na));
ans = {};
ans["K+"]  =  2*rho;
ans["Na+"] = -3*rho;
return ans;
end

-- AE
K_AE = 10e-3;
function AE(t,zeta,S,S_out)
local bicarb = S["HCO3-"];
local rho = bicarb/(bicarb+K_AE);
ans = {}
ans["HCO3-"] = -rho;
ans["Cl-"]   =  rho;
return ans;
end

-- NHE
L0    = 1000;
Kr    = 1.8e-8;
Kt    = 3.6e-6;
KNae  = 31e-3;
KHout = 2e-7;
function NHE(t,zeta,S,S_out)
local C    = Kr/Kt;
local x    = S["H+"]/Kr;
local xp1  = 1+x;
local cxp1 = 1+C*x;
local sig_num = x*xp1 + L0*C*x*cxp1;
local sig_den = L0*cxp1^2 + xp1^1;
local sig     = sig_num/sig_den;
local Nae     = S_out["Na+"];
local KNaeEff = KNae * (1+S_out["H+"]/KHout);
local sig_out = Nae / (KNaeEff+Nae);
sig = sig * sig_out;
ans = {}
ans["H+"]  = -sig;
ans["Na+"] =  sig;
return ans;
end

t_run      = 60*60;
dt         = 0.1;
save_every = 10;


