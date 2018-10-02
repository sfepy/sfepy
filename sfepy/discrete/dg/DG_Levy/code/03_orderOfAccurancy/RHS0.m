function L = RHS0(a, x1, x2, c, c_minus, c_plus, l)

int0 = 0; 

al = (x2 - x1) / (2 * l + 1);

tok_minus = a * ((c_minus(1)) + (c(1))) / 2 ...
     + abs(a) * ((c_minus(1)) - (c(1))) / 2;
      
tok_plus = a * ((c(1)) + (c_plus(1))) / 2 ...
    + abs(a) * ((c(1)) - (c_plus(1))) / 2;

Lal = a * (c(1) * int0) ...
      + tok_minus * (- 1) ^ l - tok_plus;
L = Lal / al;
end