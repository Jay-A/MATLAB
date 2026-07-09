clear all, close all

% load('L2_Errors.mat');
% B = L2_Errors;
% H = (1/2).^(B(2:end,1));
% 
% loglog(H , B(2:end,4),'-kd', 'MarkerEdgeColor', 'k', ...
%     'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);

% xlim([10^(0) 10^1]);
% ylim([10^(-12) 1]);

load('L2_Errors.mat');
L2 = L2_Errors;
H = (1/2).^L2_Errors(2:end,1);

p1 = 3; p2 = 10;
h1 = 0; h2 = 4;

p=3; loglog(H, L2(2:end, p-1),'--ks', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
hold on
p=4; loglog(H, L2(2:end, p-1),'-kd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
p=5; loglog(H, L2(2:end, p-1),'-.ko', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
p=6; loglog(H, L2(2:end, p-1),'-kp', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
p=7; loglog(H, L2(2:end, p-1),'--kh', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
p=8; loglog(H, L2(2:end, p-1),'-k^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
p=9; loglog(H, L2(2:end, p-1),'-.kx', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);
p=10; loglog(H, L2(2:end, p-1),'-k>', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.8 .8 .8], 'LineWidth',1.2, 'MarkerSize', 7);

xlim([10^(-1.95) 1.25]);
ylim([10^(-12) 40.25]);
plot([1e-3, 10], [1e-5, 1e-5], '--k')

Os = [polyfit(log(H(2:6)), log(L2_Errors(3:7,3-1)),1);
polyfit(log(H(2:6)), log(L2_Errors(3:7,4-1)),1);
polyfit(log(H(2:6)), log(L2_Errors(3:7,5-1)),1);
polyfit(log(H(2:6)), log(L2_Errors(3:7,6-1)),1);
polyfit(log(H(2:6)), log(L2_Errors(3:7,7-1)),1);
polyfit(log(H(1:5)), log(L2_Errors(2:6,8-1)),1);
polyfit(log(H(1:4)), log(L2_Errors(2:5,9-1)),1);
polyfit(log(H(1:4)), log(L2_Errors(2:5,10-1)),1)];

eval(sprintf("legend( {'p=3:  Convergence ~O(%.1f)', 'p=4:  Convergence ~O(%.1f)', 'p=5:  Convergence ~O(%.1f)', 'p=6:  Convergence ~O(%.1f)', 'p=7:  Convergence ~O(%.1f)', 'p=8:  Convergence ~O(%.1f)', 'p=9:  Convergence ~O(%.1f)', 'p=10:  Convergence ~O(%.1f)'},'FontSize',15, 'Location', 'southeast','position',[0.245 0.675 0.3 0.2]);", Os(1,1), Os(2,1), Os(3,1), Os(4,1), Os(5,1), Os(6,1), Os(7,1), Os(8,1)));

xlabel('$\log(h)$','Interpreter','latex','FontSize',24)
% ylabel('$\log\left( \left\Vert e \right\Vert_{L_2} \right)$','Interpreter','latex','FontSize',17)
ylabel('$\log\left(\sqrt{\int_\Omega\, \left(f-u\right)^2\, d\Omega}\right)$','Interpreter','latex','FontSize',24)

x0 = 10; y0 = 10;
width = 900;
height = 400;
set(gcf,'units','points','position',[x0,y0,width,height])

hold off