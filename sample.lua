
Kw = 1e-14;


species =
{
    { "H+",   1 },
    { "HO-", -1 },
    { "K+",   1 , "SP_K"}
};


equilibria =
{
    { "water", Kw, { 1, "H+" }, {1,"HO-" } }
}

function SP_K(t,zeta)
return zeta;
end
