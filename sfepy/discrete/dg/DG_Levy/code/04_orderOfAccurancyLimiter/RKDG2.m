function [W, xx, t] = RKDG2(e, a, T, int_p, int_l, pp, limiter)

% /~~~~~~~~~~~~~~~~~~ DISKRETIZACE V PROSTORU A ČASE ~~~~~~~~~~~~~~~~~~~\ %
dx = (int_p - int_l) / e;           % Krok
X = linspace(int_l,int_p,e + 1)';   % Souřadnice prvků
% CFL podmínka:
dt = dx / (abs(a) * 5); % adt/dx < 1/(2k+1), kde k je nejvyšší řád 
                        % užitých polynomů
t = zeros(ceil(T / dt) + 1,1);
fprintf('Bude provedeno %d časových kroků\n\n', ceil(T / dt));
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~ INICIALIZACE MATIC ~~~~~~~~~~~~~~~~~~~~~~~~~\ %
Wl = zeros(e,1);
Wr = zeros(e,1);
C = zeros(e,3);
C1 = zeros(e,3);
C2 = zeros(e,3);
C3 = zeros(e,3);
W = zeros(2 * e,1);
xx = zeros(2 * e,1);
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~ POČÁTEČNÍ PODMÍNKA ~~~~~~~~~~~~~~~~~~~~~~~~~\ %
Wl(:,1) = pp(X(1:end - 1));
Wr(:,1) = pp(X(2:end));
for n = 1:e
    psi0 = @(x) 1;
    psi1 = @(x) 2 * (x - (X(n + 1) + X(n)) / 2) / (X(n + 1) - X(n));
    psi2 = @(x) 6 * (x - (X(n + 1) + X(n)) / 2) ^ 2 / (X(n + 1) - X(n)) ^ 2 - 1 / 2;
    int0 = Gauss(X(n),X(n + 1) - X(n), @(x) pp(x) * psi0(x));
    int1 = Gauss(X(n),X(n + 1) - X(n), @(x) pp(x) * psi1(x));
    int2 = Gauss(X(n),X(n + 1) - X(n), @(x) pp(x) * psi2(x));
    C(n,1) = int0 * (2 * 0 + 1) / (X(n + 1) - X(n));
    C(n,2) = int1 * (2 * 1 + 1) / (X(n + 1) - X(n));
    C(n,3) = int2 * (2 * 2 + 1) / (X(n + 1) - X(n));
end
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~ RUNGE-KUTTA + DISCONTINUOUS GALERKIN ~~~~~~~~~~~~~~~~\ %
for j = 2:length(t)
   % /========================== RUNGE - KUTTA =========================\ %
   % 1. STUPEŇ:
   for n = 1:e
       if n == 1    % Ošetření okrajových podmínek (periodicity)
           m = e;     o = n + 1;
       elseif n == e
           m = n - 1; o = 1;
       else
           m = n - 1; o = n + 1;
       end
       C1(n,1) = C(n,1) + dt * RHS2(a,X(n),X(n + 1),C(n,:),C(m,:), C(o,:),0);
       C1(n,2) = C(n,2) + dt * RHS2(a,X(n),X(n + 1),C(n,:),C(m,:), C(o,:),1);
       C1(n,3) = C(n,3) + dt * RHS2(a,X(n),X(n + 1),C(n,:),C(m,:), C(o,:),2);
   end
   if strcmp(limiter, 'on') == 1
       C1 = Limiter(C1, 3);  
   end
   % 2. STUPEŇ:
   for n = 1:e
       if n == 1    % Ošetření okrajových podmínek (periodicity)
           m = e;     o = n + 1;
       elseif n == e
           m = n - 1; o = 1;
       else
           m = n - 1; o = n + 1;
       end
       C2(n,1) = 3 * C(n,1) / 4 + C1(n,1) / 4 + ... 
                 dt * RHS2(a,X(n),X(n + 1),C1(n,:),C1(m,:),C1(o,:),0) / 4;
       C2(n,2) = 3 * C(n,2) / 4 + C1(n,2) / 4 + ...
                 dt * RHS2(a,X(n),X(n + 1),C1(n,:),C1(m,:),C1(o,:),1) / 4;
       C2(n,3) = 3 * C(n,3) / 4 + C1(n,3) / 4 + ...
                 dt * RHS2(a,X(n),X(n + 1),C1(n,:),C1(m,:),C1(o,:),2) / 4;
   end
   if strcmp(limiter, 'on') == 1
       C2 = Limiter(C2, 3);  
   end
   % 3. STUPEŇ:
   for n = 1:e
       if n == 1    % Ošetření okrajových podmínek (periodicity)
           m = e;     o = n + 1;
       elseif n == e
           m = n - 1; o = 1;
       else
           m = n - 1; o = n + 1;
       end
       C3(n,1) = C(n,1) / 3 + 2 * C2(n,1) / 3 + ...
                 2 * dt * RHS2(a,X(n),X(n + 1),C2(n,:),C2(m,:),C2(o,:),0) / 3;
       C3(n,2) = C(n,2) / 3 + 2 * C2(n,2) / 3 + ...
                 2 * dt * RHS2(a,X(n),X(n + 1),C2(n,:),C2(m,:),C2(o,:),1) / 3;
       C3(n,3) = C(n,3) / 3 + 2 * C2(n,3) / 3 + ...
                 2 * dt * RHS2(a,X(n),X(n + 1),C2(n,:),C2(m,:),C2(o,:),2) / 3;
   end
   if strcmp(limiter, 'on') == 1
       C3 = Limiter(C3, 3);  
   end      
   % \==================================================================/ %

   % /========== ULOŽENÍ HODNOTY ŘEŠENÍ NA UZLU ZLEVA A ZPRAVA =========\ %
   C = C3;
   for n = 1:e
       Wl(n) = C(n,1) - C(n,2) + C(n,3);
       Wr(n) = C(n,1) + C(n,2) + C(n,3);
   end
   % \==================================================================/ %

   % /===================== ULOŽENÍ AKTUÁLNÍHO ČASU ====================\ %
   t(j) = t(j - 1) + dt;
   if ( mod(j - 1,50) ) == 0
       fprintf('%d. krok výpočtu:   t = %1.3f\n', j - 1, t(j));
   end
   % \==================================================================/ %
end
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VÝSTUPY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\ %
xx(1) = X(1);
xx(end) = X(end);
xx(2:2:end - 1) = X(2:end - 1);
xx(3:2:end - 1) = X(2:end - 1);     % Síť pro DGFEM

W(1) = Wl(1);
W(end) = Wr(end);
W(2:2:end - 1) = Wr(1:end - 1);
W(3:2:end - 1) = Wl(2:end);
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %