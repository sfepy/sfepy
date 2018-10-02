function L = RHS1(a, x1, x2, c, c_minus, c_plus, l)

switch l
    case 0
        int0 = 0; int1 = 0; 
    case 1
        int0 = 2; int1 = 0; 
end

al = (x2 - x1) / (2 * l + 1);

tok_minus = a * ((c_minus(1) + c_minus(2)) + (c(1) - c(2))) / 2 ...
     + abs(a) * ((c_minus(1) + c_minus(2)) - (c(1) - c(2))) / 2;
      
tok_plus = a * ((c(1) + c(2)) + (c_plus(1) - c_plus(2))) / 2 ...
    + abs(a) * ((c(1) + c(2)) - (c_plus(1) - c_plus(2))) / 2;

Lal = a * (c(1) * int0 + c(2) * int1) ...
      + tok_minus * (- 1) ^ l - tok_plus;
L = Lal / al;
end