clear all; clc; format long;
tic
% /~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETRY ÚLOHY ~~~~~~~~~~~~~~~~~~~~~~~~~~\ %
e = 16;            % Poèet prvkù
a = 1;            % Rychlost
T = 1;              % Maximální èas
% Pravý kraj intervalu je 1
% Levý kraj intervalu je 0
% \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/ %

% Meøení numerického øádu pøesnosti probíhá na 5 sítích, na každé je
%   dvojnásobný poèet prvkù oproti pøedchozí.
% Pro urèení øádu je brána max norma a L2 norma z rozdílu pøesného øešení 
%   a numerického øešení.

norma_max0 = zeros(5,1);
norma_max1 = zeros(5,1);
norma_max2 = zeros(5,1);
norma_L20 = zeros(5,1);
norma_L21= zeros(5,1);
norma_L22= zeros(5,1);
h = zeros(5,1);
for i = 1:5
    
    h(i) = 1/(e);
    fprintf('\n===========================================================\n');
    fprintf('            Výpoèet s %d prvky, tj. h = %.5f\n', e, h(i));
    fprintf('===========================================================\n\n');
    X = linspace(0,1,e + 1)';
    % POÈÁTEÈNÍ PODMÍNKA ... SINUS:
    pp = @(x) sin(x * 2 * pi);

    fprintf('   × k = 0 \t');
    [W_num_0,xx,t] = RKDG0(e, a, T, 1, 0, pp);
    fprintf('\n   × k = 1 \t');
    [W_num_1,~,~] = RKDG1(e, a, T, 1, 0, pp);
    fprintf('\n   × k = 2 \t');
    [W_num_2,~,~] = RKDG2(e, a, T, 1, 0, pp);

    % ANALYTICKÉ ØEŠENÍ ... SINUS:
    W_an = sin( pi .* (xx - a * t(end)) ./ (1 / 2) );
    
    
    figure(i)
    plot(xx,W_an,'y','LineWidth',3); hold on;
    plot(xx,W_num_0,'b',xx,W_num_1,'k',xx,W_num_2,'r')
    ylim([min(W_an) - 0.1, max(W_an) + 0.1]); 
    xlabel('$x$', 'Interpreter', 'latex'); 
    ylabel('$w(x,1)$',  'Interpreter', 'latex'); 
    leg1 = legend({'$\mathrm{exact}$','$k = 0$', '$k = 1$', '$k = 2$'},...
              'Location', 'northeast');
    set(leg1,'Interpreter','latex');
    set(gca, 'TickLabelInterpreter', 'latex')

    
    norma_max0(i) = max(abs(W_num_0 - W_an));
    norma_max1(i) = max(abs(W_num_1 - W_an));
    norma_max2(i) = max(abs(W_num_2 - W_an));
    
    mean_an = zeros(1,e);
    mean_num0 = zeros(1,e);
    mean_num1 = zeros(1,e);
    mean_num2 = zeros(1,e);
    for j = 1:e
        mean_an(j) = (W_an(j + 1) + W_an(j)) / 2;
        mean_num0(j) = (W_num_0(j + 1) + W_num_0(j)) / 2;
        mean_num1(j) = (W_num_1(j + 1) + W_num_1(j)) / 2;
        mean_num2(j) = (W_num_2(j + 1) + W_num_2(j)) / 2;
    end
    norma_L20(i) = (1 / e * sum(abs(mean_num0 - mean_an) .^ 2)) ^ (1 / 2);
    norma_L21(i) = (1 / e * sum(abs(mean_num1 - mean_an) .^ 2)) ^ (1 / 2);
    norma_L22(i) = (1 / e * sum(abs(mean_num2 - mean_an) .^ 2)) ^ (1 / 2);
    
    e = e * 2;
end

fprintf('\n===========================================================\n');
fprintf('                  Numerický øád pøesnosti:\n');
fprintf('===========================================================\n\n');
for i = 2:5
    p0 = log2(norma_max0(i - 1) / norma_max0(i));
    p1 = log2(norma_max1(i - 1) / norma_max1(i));
    p2 = log2(norma_max2(i - 1) / norma_max2(i));
    q0 = log2(norma_L20(i - 1) / norma_L20(i));
    q1 = log2(norma_L21(i - 1) / norma_L21(i));
    q2 = log2(norma_L22(i - 1) / norma_L22(i));
    if i == 2
        e = 16;
        fprintf('+-------+------ k = 0 -------+------ k = 1 -------+------ k = 2 -------+\n');
        fprintf('|   e   |  ||E||_L2  |   p   |  ||E||_L2  |   p   |  ||E||_L2  |   p   |\n');
        fprintf('+-------+--------------------+--------------------+--------------------+\n');
        fprintf('|  %d   | %.8f |   -   | %.8f |   -   | %.8f |   -   |\n',...
                e,norma_L20(i - 1), norma_L21(i - 1), norma_L22(i - 1));
    end
    e = e * 2;
    if length(num2str(e)) == 3
        fprintf('|  %d  | %.8f | %.3f | %.8f | %.3f | %.8f | %.3f |\n',...
                e,norma_L20(i), q0, norma_L21(i), q1,norma_L22(i), q2);
    else
        fprintf('|  %d   | %.8f | %.3f | %.8f | %.3f | %.8f | %.3f |\n',...
                e,norma_L20(i), q0, norma_L21(i), q1,norma_L22(i), q2);
    end
    if i == 5
        fprintf('+-------+--------------------+--------------------+--------------------+\n\n');
    end
end


figure(i + 1)
h1 = axes;
loglog(h,norma_L20,'-bo',...
       h,norma_L21,'-ko',...
       h,norma_L22,'-ro','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','auto');
set(h1, 'Xdir', 'reverse');
title('Chyba aproximace', 'Interpreter', 'latex');
xlabel('$\Delta x$', 'Interpreter', 'latex'); 
ylabel('$\left \| w - w_h \right \|_{L^2}$',  'Interpreter', 'latex'); 
leg1 = legend({'$k = 0$','$k = 1$','$k = 2$'},'Location', 'southwest');
set(leg1,'Interpreter','latex');
grid on; 
axis square;
set(gca, 'TickLabelInterpreter', 'latex')
toc





