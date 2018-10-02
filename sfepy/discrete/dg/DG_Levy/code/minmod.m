function M = minmod(a,b,c)

if sign(a) == sign(b) && sign(b) == sign(c)
    M = sign(a) * min([abs(a), abs(b), abs(c)]);
else
    M = 0;
end