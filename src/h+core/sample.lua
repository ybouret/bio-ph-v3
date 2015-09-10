local exp = math.exp;
-- -----------------------------------------------------------------------------
--
-- Chemistry
--
-- -----------------------------------------------------------------------------

-- -----------------------------------------------------------------------------
-- components
-- -----------------------------------------------------------------------------

lib =
{
    {"H+",   1},
    {"HO-", -1},
    {"Na+",  1},
    {"Cl-", -1},
    {"OSM",  0}
};

-- -----------------------------------------------------------------------------
-- reactions
-- -----------------------------------------------------------------------------
eqs =
{
    { "water", 1e-14, { 1, "H+"}, {  1, "HO-"} },
};



-- -----------------------------------------------------------------------------
-- inside
-- -----------------------------------------------------------------------------
ini =
{
    { "E/N" },
    { "osmolarity", 300e-3 },
    { 10e-3, {1,"Na+"} },
    { 20e-3, {1,"Cl-"} }
};

-- -----------------------------------------------------------------------------
-- outside
-- -----------------------------------------------------------------------------
out0 =
{
    { "E/N" },
    { "osmolarity", 300e-3},
    { 140e-3, {1,"Na+"} },
    { 100e-3, {1,"Cl-"} }
};

out =
{
    "out0"
};

function weights(t)
return 1;
end

-- -----------------------------------------------------------------------------
--
-- Physics
--
-- -----------------------------------------------------------------------------

T       = 298;
Cm      = 1;   -- pF/cm^2

Capa_exp = 207e-12;               -- from article 207 pF
Surf_exp = Capa_exp / (1e-14*Cm); -- in micron^2

-- -----------------------------------------------------------------------------
--
-- effectors
--
-- -----------------------------------------------------------------------------
eff =
{
     "lambda_Na" -- ,"lambda_Cl"
    --"lambda_Cl","lambda_Na"
};

-- -----------------------------------------------------------------------------
--
-- parameters
--
-- -----------------------------------------------------------------------------
Em    = 0;
zeta0 = Em * F/(R*T);

a    = 12.5; -- in microns
b    = 5;    -- in microns
c    = 5;    -- in microns

surface = EllipsoidSurface(a,b,c);
volume  = EllipsoidVolume(a,b,c);


-- -----------------------------------------------------------------------------
--
-- solver and run
--
-- -----------------------------------------------------------------------------
ftol   = 1e-5;
diff_h = 1e-6;


dt      = 0.1;
dt_save = 0.5;
t_run   = 1000;

-- -----------------------------------------------------------------------------
-- Surface * Permeability in micron^3 / s, fitted for Surf_exp
-- -----------------------------------------------------------------------------

function SP_K( t, x )
return ((3.70714279) * exp( (0.4310085034) * x ) + (1.507274188));
end

-- -----------------------------------------------------------------------------
-- Surface * Permeability in micron^3 / s, fitted for Surf_exp
-- -----------------------------------------------------------------------------
function SP_Cl(t,x)
return (1.194531961) * exp( (0.05658342092) * x );
end

-- -----------------------------------------------------------------------------
-- Surface * Permeability in micron^3 / s, fitted for Surf_exp
-- -----------------------------------------------------------------------------

function SP_Na(t,x)
return (1.360058333) * exp( (0.05223441992) * x );
end


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
local Perm = SP_Cl(t,zeta)/Surf_exp; -- in microns/s
local Flux = Perm * Psi(zz)*(Cout[Cl]-Cin[Cl]*exp(zz)); -- in moles/L*microns/s
local J    = 1e-3*Flux;                                 -- in moles/m^2/s
a[Cl]      = J;
return a;

end
