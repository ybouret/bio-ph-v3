-- -----------------------------------------------------------------------------
--
-- times
--
-- -----------------------------------------------------------------------------
ftol   = 1e-4; -- differential fractional tolerance
diff_h = 1e-5; -- initial adaptive time step between to time steps

dt      = 0.1;
dt_save = 0.2;
t_run   = 15*60;

lib =
{
    {"H+",1},
    {"HO-",-1},
    {"OSM",0}
};

eqs =
{
    { "water", 1e-14, { 1, "H+"}, {  1, "HO-"} },
};

eff =
{

};



ini =
{
    { "E/N" },
    { "osmolarity", 300e-3 }
};

out0 =
{
    { "E/N" },
    { "osmolarity", 300e-3}
};

out =
{
    "out0"
};

function weights(t)
return 1;
end

volume =1;
surface=1;

-- -----------------------------------------------------------------------------
-- Common surface capacitance
-- -----------------------------------------------------------------------------


Cm = 10 * 1.0e-15; -- Farad/ micron^2
