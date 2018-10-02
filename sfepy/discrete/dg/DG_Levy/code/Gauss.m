function int = Gauss(a,h,f)
%--------------------------------------------------------------------------
%   Dvobodový Gaussùv kvadraturní vzorec transformovaný 
%   z intervalu <-1,1> na <a,b>
%       a ... levá mez intervalu <a,b>
%       h ... délka intervalu (b-a)
%       f ... integrovaná funkce
%--------------------------------------------------------------------------

x_1 = a + h*(-sqrt(1/3)+1)/2;
x_2 = a + h*(sqrt(1/3)+1)/2;        %body
w_1 = h/2;
w_2 = h/2;                          %váhy

int = w_1*f(x_1) + w_2*f(x_2);      %hodnota integrálu

end