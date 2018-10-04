function [C] = Limiter(c, l)

e = size(c,1);
C = c;

for i = 1:e
    if i == 1    % Ošetření okrajových podmínek (periodicity)
        m = e;     o = i + 1;
    elseif i == e
        m = i - 1; o = 1;
    else
        m = i - 1; o = i + 1;
    end
    for k = l:-1:2
        
        C(i,k) = minmod(c(i,k), c(o, k - 1) - c(i, k - 1), c(i, k - 1) - c(m, k - 1));

        if C(i,k) == c(i,k) 
            break
        end
    end
end


