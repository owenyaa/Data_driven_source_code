%% Author : TAO ZHANG  * zt1996nic@gmail.com *
% Created Time : 2023-07-20 00:00
% Last Revised : TAO ZHANG ,2023-08-01 00:00
% Remark : Fractional Chaotic System: fractional order Chua’s circuit

clear;clc;close all;
addpath('./utils');

% Fractional order
q1=0.98; 
% Time setting
h=0.001;t0=0;tfinal=100;
% Initial condition
y0=[0.1; 0.1; 0.1]; 
% y0=[0.2; -0.1; 0.1];
% Add noise
fr = 0.01;
% The order setting of the library
X_OrderMax=2;
Trig_OrderMax=0;
nonsmooth_OrderMax=0;
% Algorithm parameter settings
params_alg.maxit = 200; 
params_alg.del = 20; 
params_alg.dq = 0.001;
params_alg.q0 = 0.9;
params_alg.lambda = 0.05;

% FDE Slover
[t, y] = fde12(q1,@(t, y)SimpleChua(t, y),t0,tfinal,y0,h);
x=y';
[N, DIM] = size(x);


% % True data
% figure(1);
% plot3(x(:,1),x(:,2),x(:,3),'r-', 'LineWidth', 1.5);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
% view(27,16)
% hold off;

% % Add noise
xn = x + fr * (2 * rand(N, DIM) - 1);
X=xn;

tic;
% Construct Library
[Theta,Sym]=LIBB(X,X_OrderMax,Trig_OrderMax,nonsmooth_OrderMax);
% [Theta,Sym]=LIBC(X,X_OrderMax,Trig_OrderMax,nonsmooth_OrderMax);
% Iterative hard thresholding
[Xi,q,outputs,options] = IHT_FO(X, h,t0, tfinal, Theta, params_alg);


% Minimized selection By model err
distancesM = sqrt(((min(outputs.ObFu1) - outputs.ObFu1)).^2 + ...
    ((min(outputs.ObFu3) - outputs.ObFu3)).^2);
[min_distance_M, min_index_M] = min(distancesM);
q_optM = outputs.q(min_index_M);
Xi_optM=outputs.Xi{min_index_M};
RXi_optM = Xi_optM(:);
toc;
TF=toc;
% Recovery model By model err
[~, SM] = fde12(q_optM,@(t, y)sparse_model(t, y, RXi_optM),t0,tfinal,x(1,:)',h);
% Output result
disp("minimum distance by model error: " + min_distance_M);
disp("The serial number of the corresponding q: " + min_index_M);
disp("The corresponding q value: " + q_optM);
disp("Dictionary Index: " + Sym' + "     |      Coefficient Value " + num2str(Xi_optM(1:end,:)));


% Minimized selection By response err
distancesX = sqrt(((min(outputs.ObFu1) - outputs.ObFu1)).^2 + ...
    ((min(outputs.ObFu2) - outputs.ObFu2)).^2);
[min_distance_X, min_index_X] = min(distancesX);
q_optX = outputs.q(min_index_X);
Xi_optX=outputs.Xi{min_index_X};
RXi_optX = Xi_optX(:);
% Recovery model By response err
[~, SX] = fde12(q_optX,@(t, y)sparse_model(t, y, RXi_optX),t0,tfinal,x(1,:)',h);
% Output result
disp("minimum distance by response error: " + min_distance_X);
disp("The serial number of the corresponding q: " + min_index_X);
disp("The corresponding q value: " + q_optX);
disp("Dictionary Index: " + Sym' + "     |      Coefficient Value " + num2str(Xi_optX(1:end,:)));

% Minimized selection By OB2 OB3
distancesS = sqrt(((min(outputs.ObFu2) - outputs.ObFu2)).^2 + ...
    ((min(outputs.ObFu3) - outputs.ObFu3)).^2);
[min_distance_S, min_index_S] = min(distancesS);
q_optS = outputs.q(min_index_S);
Xi_optS=outputs.Xi{min_index_S};
RXi_optS = Xi_optS(:);
% Recovery model By response err
[~, SS] = fde12(q_optS,@(t, y)sparse_model(t, y, RXi_optS),t0,tfinal,x(1,:)',h);
% Output result
disp("minimum distance by S: " + min_distance_S);
disp("The serial number of the corresponding q: " + min_index_S);
disp("The corresponding q value: " + q_optS);
disp("Dictionary Index: " + Sym' + "     |      Coefficient Value " + num2str(Xi_optS(1:end,:)));


% Minimized selection By OB1 OB2 OB3
distancesY = sqrt(((min(outputs.ObFu1) - outputs.ObFu1)).^2 + ...
    ((min(outputs.ObFu2) - outputs.ObFu2)).^2 + ...
    ((min(outputs.ObFu3) - outputs.ObFu3)).^2);
[min_distance_Y, min_index_Y] = min(distancesY);
q_optY = outputs.q(min_index_Y);
Xi_optY=outputs.Xi{min_index_Y};
RXi_optY = Xi_optY(:);
% Recovery model By response err
[~, SY] = fde12(q_optY,@(t, y)sparse_model(t, y, RXi_optY),t0,tfinal,x(1,:)',h);
% Output result
disp("minimum distance by Y: " + min_distance_Y);
disp("The serial number of the corresponding q: " + min_index_Y);
disp("The corresponding q value: " + q_optY);
disp("Dictionary Index: " + Sym' + "     |      Coefficient Value " + num2str(Xi_optY(1:end,:)));

% % True Xi

% Xi_true=[2.95	0	0
%             0    0    0
%             0    0    0
%             -2.95	0	0
%             0    0    0
%             0    0    0            
%             -3.2	1	0
%             10	-1	-14.87
%             0	1	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0];

% Xi_true=[2.95	0	0
%             -2.95	0	0
%             -3.2	1	0
%             10	-1	-14.87
%             0	1	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0];

Xi_true=[2.95	0	0
            0  0  0
            0  0  0
            -2.95	0	0
            0  0  0
            0  0  0
            -3.2	1	0
            10	-1	-14.87
            0	1	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0];

RE = (Xi_optM(abs(Xi_true)>0) - Xi_true(abs(Xi_true)>0))./Xi_true(abs(Xi_true)>0);
MRE = max(abs(RE))*100;
fprintf('The maximum relative error %2.2f\n',MRE);
% 
REX = (Xi_optX(abs(Xi_true)>0) - Xi_true(abs(Xi_true)>0))./Xi_true(abs(Xi_true)>0);
MREX = max(abs(REX))*100;
fprintf('The maximum relative error %2.2f\n',MREX);

% True data
h(1)=figure;
plot3(x(:,1),x(:,2),x(:,3),'r-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
view(27,16)
hold off;

% figure(11);
% plot3(x(50000:end,1),x(50000:end,2),x(50000:end,3),'r-', 'LineWidth', 1.5);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
% view(27,16)
% hold off;

% Measured data
h(2)=figure;
plot3(X(:,1),X(:,2),X(:,3),'k-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
view(27,16)
hold off;

% Recovery data
h(3)=figure;
plot3(x(:,1),x(:,2),x(:,3),'r-', 'LineWidth', 1.5);
hold on;
plot3(SX(1,:),SX(2,:),SX(3,:),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
h1=legend('Exact','SIFCS','FontSize',15);
set(h1,'FontName', 'Times New Roman');
title('$err_x$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
view(27,16)
hold off;

h(4)=figure;
plot3(x(:,1),x(:,2),x(:,3),'r-', 'LineWidth', 1.5);
hold on;
plot3(SM(1,:),SM(2,:),SM(3,:),'b-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
h1=legend('Exact','SIFCS','FontSize',15);
set(h1,'FontName', 'Times New Roman');
title('$err_m$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
view(27,16)
hold off;

h(41)=figure;
plot3(X(:,1),X(:,2),X(:,3),'k-', 'LineWidth', 1.5);
hold on;
plot3(SM(1,:),SM(2,:),SM(3,:),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
h1=legend('Measured','SIFCS','FontSize',15);
set(h1,'FontName', 'Times New Roman');
title('$err_m$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
view(27,16)
hold off;

h(5)=figure;
subplot(3, 1, 1);
plot(t,x(1:end,1),'r-', 'LineWidth', 1.5);
hold on;
plot(t,SX(1,1:end),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
title('$err_x$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
hold on;
% axis tight; 
subplot(3, 1, 2);
plot(t,x(1:end,2),'r-', 'LineWidth', 1.5);
hold on;
plot(t,SX(2,1:end),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% axis tight; 
subplot(3, 1, 3);
plot(t,x(1:end,3),'r-', 'LineWidth', 1.5);
hold on;
plot(t,SX(3,1:end),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
% axis tight; 
hold off;

% % YY=0.5*abs(x(1:end,2)+1)-0.5*abs(x(1:end,2)-1);
% % XX=0.5*abs(x(1:end,1)+1)-0.5*abs(x(1:end,1)-1);
% % ZZ=0.5*abs(x(1:end,3)+1)-0.5*abs(x(1:end,3)-1);
% % h(51)=figure;
% % subplot(3, 1, 1);
% % plot(t,x(1:end,1),'r-', 'LineWidth', 1.5);
% % hold on;
% % plot(t,XX,'k--', 'LineWidth', 1.5);
% % set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% % h1=legend('$x$','$0.5|x+1|-0.5|x-1|$','FontSize',15);
% % set(h1, 'Interpreter', 'latex');
% % ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% % hold on;
% % % axis tight; 
% % subplot(3, 1, 2);
% % plot(t,x(1:end,2),'r-', 'LineWidth', 1.5);
% % hold on;
% % plot(t,YY,'k--', 'LineWidth', 1.5);
% % set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% % h1=legend('$y$','$0.5|y+1|-0.5|y-1|$','FontSize',15);
% % set(h1, 'Interpreter', 'latex');
% % ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% % % yline(0,'k--','LineWidth',2.5);
% % % axis tight; 
% % subplot(3, 1, 3);
% % plot(t,x(1:end,3),'r-', 'LineWidth', 1.5);
% % hold on;
% % plot(t,ZZ,'k--', 'LineWidth', 1.5);
% % set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% % h1=legend('$z$','$0.5|z+1|-0.5|z-1|$','FontSize',15);
% % set(h1, 'Interpreter', 'latex');
% % ylabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
% % % axis tight; 
% % hold off;

% % YY=0.5*abs(x(1:end,2)+1)-0.5*abs(x(1:end,2)-1);
% % h(52)=figure;
% % plot(t,YY,'r-', 'LineWidth', 1.5);
% % hold on;
% % plot(t,x(1:end,2),'k--', 'LineWidth', 1.5);
% % set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% % ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% % % yline(0,'k--','LineWidth',2.5);


h(6)=figure;
subplot(3, 1, 1);
plot(t,x(1:end,1),'r-', 'LineWidth', 1.5);
hold on;
plot(t,SM(1,1:end),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
title('$err_m$', 'Interpreter', 'latex', 'FontSize', 20);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
hold on;
% axis tight; 
subplot(3, 1, 2);
plot(t,x(1:end,2),'r-', 'LineWidth', 1.5);
hold on;
plot(t,SM(2,1:end),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% axis tight; 
subplot(3, 1, 3);
plot(t,x(1:end,3),'r-', 'LineWidth', 1.5);
hold on;
plot(t,SM(3,1:end),'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
% axis tight; 
hold off;

% q Selection graph 3D: Color
h(7)=figure;
colormap('jet');
scatter3(outputs.q, outputs.ObFu1, outputs.ObFu2, 80, distancesX, 'filled');
colorbar;
hold on;
scatter3(outputs.q(min_index_X), outputs.ObFu1(min_index_X), outputs.ObFu2(min_index_X), ...
    300, 'p', 'filled', 'MarkerEdgeColor',...
    'k', 'MarkerFaceColor', 'm');
x_val = outputs.q(min_index_X);
y_val = outputs.ObFu1(min_index_X);
z_val = outputs.ObFu2(min_index_X);
plot3([x_val, x_val], [y_val, y_val], [0, z_val], '--k', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$q$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| \Xi \|_0$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$err_x$', 'Interpreter', 'latex', 'FontSize', 18);
hold off;


h(8)=figure;
colormap('jet');
scatter3(outputs.q, outputs.ObFu1, outputs.ObFu3, 80, distancesM, 'filled');
colorbar;
hold on;
scatter3(outputs.q(min_index_M), outputs.ObFu1(min_index_M), outputs.ObFu3(min_index_M), ...
    300, 'p', 'filled', 'MarkerEdgeColor',...
    'k', 'MarkerFaceColor', 'm');
x_val = outputs.q(min_index_M);
y_val = outputs.ObFu1(min_index_M);
z_val = outputs.ObFu3(min_index_M);
plot3([x_val, x_val], [y_val, y_val], [0, z_val], '--k', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$q$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\| \Xi \|_0$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$err_m$', 'Interpreter', 'latex', 'FontSize', 18);
hold off;

h(9)=figure;
colormap('jet');
scatter3(outputs.q, outputs.ObFu2, outputs.ObFu3, 80, distancesS, 'filled');
colorbar;
hold on;
scatter3(outputs.q(min_index_S), outputs.ObFu2(min_index_S), outputs.ObFu3(min_index_S), ...
    300, 'p', 'filled', 'MarkerEdgeColor',...
    'k', 'MarkerFaceColor', 'm');
x_val = outputs.q(min_index_S);
y_val = outputs.ObFu2(min_index_S);
z_val = outputs.ObFu3(min_index_S);
plot3([x_val, x_val], [y_val, y_val], [0, z_val], '--k', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$q$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$err_x$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$err_m$', 'Interpreter', 'latex', 'FontSize', 18);
hold off;

h(10)=figure;
colormap('jet');
scatter3(outputs.ObFu1, outputs.ObFu2, outputs.ObFu3, 80, distancesY, 'filled');
colorbar;
hold on;
scatter3(outputs.ObFu1(min_index_Y), outputs.ObFu2(min_index_Y), outputs.ObFu3(min_index_Y), ...
    300, 'p', 'filled', 'MarkerEdgeColor',...
    'k', 'MarkerFaceColor', 'm');
x_val = outputs.ObFu1(min_index_Y);
y_val = outputs.ObFu2(min_index_Y);
z_val = outputs.ObFu3(min_index_Y);
plot3([x_val, x_val], [y_val, y_val], [0, z_val], '--k', 'LineWidth', 2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
xlabel('$\| \Xi \|_0$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$err_x$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$err_m$', 'Interpreter', 'latex', 'FontSize', 18);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%