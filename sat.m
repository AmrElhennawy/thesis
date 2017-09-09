function y = sat(s)
    k = 0.15;
    y = (2./(1+exp(-k*s))) - 1.0;
end
