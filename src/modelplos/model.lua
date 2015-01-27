
lib =
{
    {"H+",  1},
    {"HO-",-1}
};

eqs =
{
    { "water", 1e-14, { 1, "H+"}, { 1, "HO-"} }
};

eff =
{
};

ini =
{
    { "E/N" }
};


out1 =
{
    { "E/N" }
};

out =
{
    "out1"
};

function weight(t)
return 1
end
