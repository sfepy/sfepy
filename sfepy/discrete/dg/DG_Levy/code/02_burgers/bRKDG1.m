function [W, xx, t] = bRKDG1(e, T, int_l, int_p, pp, limiter)

% /~~~~~~~~~~~~~~~~~~ DISKRETIZACE V PROSTORU A ÈASE ~~~~~~~~~~~~~~~~~~~\ %
dx = (int_p - int_l) / e;           % Krok
X = linspace(int_l,int_p,e + 1)';   % Souøadnice prvkù
Wmax = max(pp(X));
% CFL podmínka:
dt = 0.5 * dx / (abs(Wmax) * 5); % adt/dx < 1/(2k+1), kde k je nejvyšší øád 
                        % užitých polynomù
t = zeros(ceil(T / dt) + 1,1);
fprintf('Bude provedeno %d èasových krokù\n\n', ceil(T / dt));
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~ INICIALIZACE MATIC ~~~~~~~~~~~~~~~~~~~~~~~~~\ %
Wl = zeros(e,length(t));
Wr = zeros(e,length(t));
C = zeros(e,2);
C1 = zeros(e,2);
C2 = zeros(e,2);
C3 = zeros(e,2);
W = zeros(2 * e,length(t));
xx = zeros(2 * e,1);
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~ POÈÁTEÈNÍ PODMÍNKA ~~~~~~~~~~~~~~~~~~~~~~~~~\ %
Wl(:,1) = pp(X(1:end - 1));
Wr(:,1) = pp(X(2:end));
for n = 1:e
    psi0 = @(x) 1;
    psi1 = @(x) 2 * (x - (X(n + 1) + X(n)) / 2) / (X(n + 1) - X(n));
    int0 = Gauss(X(n),X(n + 1) - X(n), @(x) pp(x) * psi0(x));
    int1 = Gauss(X(n),X(n + 1) - X(n), @(x) pp(x) * psi1(x));
    C(n,1) = int0 * (2 * 0 + 1) / (X(n + 1) - X(n));
    C(n,2) = int1 * (2 * 1 + 1) / (X(n + 1) - X(n));
end
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~ RUNGE-KUTTA + DISCONTINUOUS GALERKIN ~~~~~~~~~~~~~~~~\ %
for j = 2:length(t)
   % /========================== RUNGE - KUTTA =========================\ %
   % 1. STUPEÒ:
   for n = 1:e
       if n == 1    % Ošetøení okrajových podmínek (periodicity)
           m = e;     o = n + 1;
       elseif n == e
           m = n - 1; o = 1;
       else
           m = n - 1; o = n + 1;
       end
       C1(n,1) = C(n,1) + dt * RHS1(Wmax,X(n),X(n + 1),C(n,:),C(m,:),C(o,:),0);
       C1(n,2) = C(n,2) + dt * RHS1(Wmax,X(n),X(n + 1),C(n,:),C(m,:),C(o,:),1);
   end
   if strcmp(limiter, 'on') == 1
       C1 = Limiter(C1, 2);  
   end
   % 2. STUPEÒ:
   for n = 1:e
       if n == 1    % Ošetøení okrajových podmínek (periodicity)
           m = e;     o = n + 1;
       elseif n == e
           m = n - 1; o = 1;
       else
           m = n - 1; o = n + 1;
       end
       C2(n,1) = 3 * C(n,1) / 4 + C1(n,1) / 4 + ... 
                 dt * RHS1(Wmax,X(n),X(n + 1),C1(n,:),C1(m,:),C1(o,:),0) / 4;
       C2(n,2) = 3 * C(n,2) / 4 + C1(n,2) / 4 + ...
                 dt * RHS1(Wmax,X(n),X(n + 1),C1(n,:),C1(m,:),C1(o,:),1) / 4;
   end
   if strcmp(limiter, 'on') == 1
       C2 = Limiter(C2, 2);  
   end
   % 3. STUPEÒ:
   for n = 1:e
       if n == 1    % Ošetøení okrajových podmínek (periodicity)
           m = e;     o = n + 1;
       elseif n == e
           m = n - 1; o = 1;
       else
           m = n - 1; o = n + 1;
       end
       C3(n,1) = C(n,1) / 3 + 2 * C2(n,1) / 3 + ...
                 2 * dt * RHS1(Wmax,X(n),X(n + 1),C2(n,:),C2(m,:),C2(o,:),0) / 3;
       C3(n,2) = C(n,2) / 3 + 2 * C2(n,2) / 3 + ...
                 2 * dt * RHS1(Wmax,X(n),X(n + 1),C2(n,:),C2(m,:),C2(o,:),1) / 3;
   end
   if strcmp(limiter, 'on') == 1
       C3 = Limiter(C3, 2);  
   end
   % \==================================================================/ %

   % /========== ULOŽENÍ HODNOTY ØEŠENÍ NA UZLU ZLEVA A ZPRAVA =========\ %
   C = C3;
   for n = 1:e
       Wl(n,j) = C(n,1) - C(n,2);
       Wr(n,j) = C(n,1) + C(n,2);
   end
   % \==================================================================/ %

   % /===================== ULOŽENÍ AKTUÁLNÍHO ÈASU ====================\ %
   t(j) = t(j - 1) + dt;
%    fprintf('%d   t = %1.3f\n', j - 1, t(j));
   % \==================================================================/ %
end
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VÝSTUPY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\ %
xx(1) = X(1);
xx(end) = X(end);
xx(2:2:end - 1) = X(2:end - 1);
xx(3:2:end - 1) = X(2:end - 1);     % Sí pro DGFEM

for j = 1:length(t)
    W(1,j) = Wl(1,j);
    W(end,j) = Wr(end,j);
    W(2:2:end - 1,j) = Wr(1:end - 1,j);
    
    W(3:2:end - 1,j) = Wl(2:end,j);
end
end


function L = RHS1(Wmax, x1, x2, c, c_minus, c_plus, l)

psi0 = @(x) 1;
psi1 = @(x) 2 * (x - (x2 + x1) / 2) / (x2 - x1);

switch l
    case 0
        psi_x = @(x) 0;
    case 1
        psi_x = @(x) 2 / (x2 - x1);
end
int = Gauss(x1, (x2 - x1), ...
      @(x) ((c(1) * psi0(x) + c(2) * psi1(x)) .^ 2) .* psi_x(x));
  
al = (x2 - x1) / (2 * l + 1);

tok_minus = 0.5 * ( 0.5 * (c_minus(1) + c_minus(2)) ^ 2 + 0.5 * (c(1) - c(2)) ^ 2 ) + ...
            0.5 * Wmax * ( (c_minus(1) + c_minus(2)) - (c(1) - c(2)) );
      
tok_plus =  0.5 * ( 0.5 * (c(1) + c(2)) ^ 2 + 0.5 * (c_plus(1) - c_plus(2)) ^ 2 ) + ...
            0.5 * Wmax * ( (c(1) + c(2)) - (c_plus(1) - c_plus(2)) );


Lal = 0.5 * int + tok_minus * (- 1) ^ l - tok_plus;
L = Lal / al;
end