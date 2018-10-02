clear all; clc;
tic
% /~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETRY ÚLOHY ~~~~~~~~~~~~~~~~~~~~~~~~~~\ %
a = 1;            % Rychlost
T = 1;              % Maximální èas
int_p = 1;          % Pravý kraj intervalu
int_l = 0;          % Levý kraj intervalu
e = 100;            % Poèet prvkù
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% /~~~~~~~~~~~~~~~~~~~~~~~~ POÈÁTEÈNÍ PODMÍNKA ~~~~~~~~~~~~~~~~~~~~~~~~~\ %
%   1) GAUSSOVA FUNKCE:
% pp = @(x) exp(-((x - (int_l + int_p) / 2).* 5 ./ ((int_l + int_p) / 2) ).^2 );

%   2) PARABOLA:
% pp = @(x) - 100 * (x - 0.1) .* (x - 0.3) .* (x >= 0.1 & x <= 0.3) + 0 .* (x < 0.1 & x > 0.3);

%   3) SINUS:
% pp = @(x) sin( x * 2 * pi / (int_p + abs(int_l)) );

%   4) JEDNOTKOVÝ SKOK:
% pp = @(x) 1 * (x >= 0.1 & x <= 0.3) + 0 * (x < 0.1 & x > 0.3);

%   5) MIX FUNKCÍ
delta = 0.005;
alpha = 25;
b = 0.75;
z = 0.15;
beta = log(4) / (36 * delta ^ 2);

ex = @(x,y) exp(- beta .* (x - y) .^ 2);
F = @(x,y) sqrt(max(1 - alpha ^ 2 * (x - y) .^ 2, 0));

pp = @(x) 0 * (x >= 0 & x < 0.1) + ...
          0.5 / 6 * (ex(x,z - delta) + ex(x,z + delta) + 4 .* ex(x,z)) .* (x >= 0.1 & x < 0.2) + ...
          0 * (x >= 0.2 & x < 0.3) + ... 
          0.5 * (x >= 0.3 & x < 0.4) + ...
          0 * (x >= 0.4 & x < 0.5) + ...
          (0.5 - abs(10 * (x - 0.55))) .* (x >= 0.5 & x < 0.6) + ...
          0 * (x >= 0.6 & x < 0.7) + ...
          0.5 / 6 * (F(x,b - delta) + F(x,b + delta) + 4 .* F(x,b)) .* (x >= 0.7 & x < 0.8) + ...
          0 * (x >= 0.8 & x < 0.9) + ...
          0 * (x >= 0.9 & x < 1);    

% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

[W, xx, t] = RKDG2(e, T, a, int_l, int_p, pp,'off');

[Wlim, ~, ~] = RKDG2(e, T, a, int_l, int_p, pp,'on');

% /~~~~~~~~~~~~~~~~~~~~~~~~~ GRAFICKÉ VÝSTUPY ~~~~~~~~~~~~~~~~~~~~~~~~~~\ %
% subplot(2,1,1)
    plot(xx,W(:,1),'-k',xx,W(:,end),'-b',xx,Wlim(:,end),'-r');
    legend({'exact','k = 2','k = 2 + limiter'},'Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex'); 
    ylabel('$w(x,1)$', 'Interpreter', 'latex'); 
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex')
% subplot(2,1,2)
%     surf(xx,t,W', 'LineStyle', 'none')
%     xlabel('$x$', 'Interpreter', 'latex'); 
%     ylabel('$t$', 'Interpreter', 'latex'); 
%     zlabel('$w(x,t)$', 'Interpreter', 'latex')
%     ylim([t(1),t(end)]);
%     xlim([xx(1),xx(end)]);
%     set(gca, 'TickLabelInterpreter', 'latex')
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %


toc





