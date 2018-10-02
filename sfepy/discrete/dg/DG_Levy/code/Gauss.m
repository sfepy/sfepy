function int = Gauss(a,h,f)
%--------------------------------------------------------------------------
%   Dvobodov� Gauss�v kvadraturn� vzorec transformovan� 
%   z intervalu <-1,1> na <a,b>
%       a ... lev� mez intervalu <a,b>
%       h ... d�lka intervalu (b-a)
%       f ... integrovan� funkce
%--------------------------------------------------------------------------

x_1 = a + h*(-sqrt(1/3)+1)/2;
x_2 = a + h*(sqrt(1/3)+1)/2;        %body
w_1 = h/2;
w_2 = h/2;                          %v�hy

int = w_1*f(x_1) + w_2*f(x_2);      %hodnota integr�lu

end