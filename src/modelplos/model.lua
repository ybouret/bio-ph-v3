
-- standard
exp  = math.exp;
tanh = math.tanh;

lib =
{
    {"H+",   1},
    {"HO-", -1},
    {"Na+",  1},
    {"Cl-", -1},
    {"K+",   1},
    {"OSM",  0}
};

eqs =
{
    { "water", 1e-14, { 1, "H+"}, { 1, "HO-"} }
};

eff =
{
    "j_Na",
    "j_K",
    "j_Cl",
    "NaK"
};


ini =
{
    { 10^(-7.2), { 1, "H+" } },
    { 140e-3,    { 1, "K+"  } },
    { 10e-3,     { 1, "Na+" } },
    { 20e-3,     { 1, "Cl-" } },
    { "osmolarity", 300e-3}
};




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

function TableSumValues(X)
ans = 0;
for key,value in pairs(X) do ans = ans+value; end
return ans;
end

function j_Na(t,Cin,Cout,params)
local zz = params["zeta"];
local Na = "Na+";

-- for key,value in pairs(Cin) do print(key,value) end
-- print(TableSumValues(Cin));

a = {};
local rho = -Psi(zz)*(Cout[Na]-Cin[Na]*exp(zz));
a[Na] = rho;
return a;
end


function j_K(t,Cin,Cout,params)
local zz = params["zeta"];
local K  = "K+";
a = {};
local rho = -Psi(zz)*(Cout[K]-Cin[K]*exp(zz));
a[K] = rho;
return a;
end

function j_Cl(t,Cin,Cout,params)
local zz = -params["zeta"];
local Cl = "Cl-";
a = {};
local rho = -Psi(zz)*(Cout[Cl]-Cin[Cl]*exp(zz));
a[Cl] = rho;
return a;
end


K_NaK = 12e-3;

function NaK(t,Cin,Cout,params)
local zeta  = params["zeta"];
local CNa   = Cin["Na+"];
local sig   = (CNa/(K_NaK+CNa)) * (1+tanh(0.39*zeta+1.28))*0.5;
a = {}
a["K+"]  =  2*sig;
a["Na+"] = -3*sig;
return a;
end



