%% Author : TAO ZHANG  * zt1996nic@gmail.com *
% Created Time : 2023-07-20 00:00
% Last Revised : TAO ZHANG ,2023-08-01 00:00
% Remark : Fractional Chaotic System: Fractional-order Simple Lorenz

clear;clc;close all;
addpath('./utils');

% Set the desired path for saving the files
% savePath = 'D:\博三\分数阶混沌_第二篇\Chinese Journal of Physics\CODE_DDPC\SimpleLorenz\c5qvar\order3\0.001noise\dq0.001';

% Fractional order
q1=0.95; 
% Syetem parameter
c=5;
% Time setting
h=0.01;t0=0;tfinal=20;
% Initial condition
y0=[1; 1; 1]; 
% y0=[-8; 7; 27];
% y0=[1; 1; 8];

% Add noise
fr = 0.1;
% The order setting of the library
X_OrderMax=2;
Trig_OrderMax=0;
nonsmooth_OrderMax=0;
% Algorithm parameter settings
params_alg.maxit = 2000; 
params_alg.del = 20; 
params_alg.dq = 0.0001;
params_alg.q0 = 0.9;
params_alg.lambda = 0.1;

% FDE Slover
[t, y] = fde12(q1,@(t, y)SimpleLorenz(t, y, c),t0,tfinal,y0,h);
x=y';
[N, DIM] = size(x);

% h(111)=figure;
% plot3(x(:,1),x(:,2),x(:,3),'r-', 'LineWidth', 1.5);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
% zlabel('$z$', 'Interpreter', 'latex', 'FontSize', 20);
% view(27,16)
% hold off;

% % Add noise
xn = x + fr * (2 * rand(N, DIM) - 1);
% xn = x + fr * randn(N, DIM); 
X=xn;

% % D^q x, in the initial value section
% DDX=10*(x(:,2)-x(:,1));
% DDY=(24-4*5)*x(:,1) - x(:,1).*x(:,3) + 5*x(:,2);
% DDZ=x(:,1).*x(:,2) - 8/3*x(:,3);
% D_qx=FO_XResponse(q1,h,x,N);
% Dq = ban(q1,N,h);
% for i=1:3
%     D_qx1(:,i) = Dq * x(:,i);
% end
% del=20;
% figure;plot(t(1:del),DDX(1:del),'r', 'LineWidth', 1.5);hold on;plot(t(1:del),D_qx(1:del,1),'k.', 'LineWidth', 1.5);
% figure;plot(t(1:del),DDY(1:del),'r', 'LineWidth', 1.5);hold on;plot(t(1:del),D_qx(1:del,2),'k.', 'LineWidth', 1.5);
% figure;plot(t(1:del),DDZ(1:del),'r', 'LineWidth', 1.5);hold on;plot(t(1:del),D_qx(1:del,3),'k.', 'LineWidth', 1.5);
% figure;plot(t(1:del),DDX(1:del),'r', 'LineWidth', 1.5);hold on;plot(t(1:del),D_qx1(1:del,1),'k.', 'LineWidth', 1.5);
% figure;plot(t(1:del),DDY(1:del),'r', 'LineWidth', 1.5);hold on;plot(t(1:del),D_qx1(1:del,2),'k.', 'LineWidth', 1.5);
% figure;plot(t(1:del),DDZ(1:del),'r', 'LineWidth', 1.5);hold on;plot(t(1:del),D_qx1(1:del,3),'k.', 'LineWidth', 1.5);
% 
% figure;plot(t(end-del:end),DDX(end-del:end),'r', 'LineWidth', 1.5);hold on;plot(t(end-del:end),D_qx(end-del:end,1),'k.', 'LineWidth', 1.5);
% figure;plot(t(end-del:end),DDY(end-del:end),'r', 'LineWidth', 1.5);hold on;plot(t(end-del:end),D_qx(end-del:end,2),'k.', 'LineWidth', 1.5);
% figure;plot(t(end-del:end),DDZ(end-del:end),'r', 'LineWidth', 1.5);hold on;plot(t(end-del:end),D_qx(end-del:end,3),'k.', 'LineWidth', 1.5);
% figure;plot(t(end-del:end),DDX(end-del:end),'r', 'LineWidth', 1.5);hold on;plot(t(end-del:end),D_qx1(end-del:end,1),'k.', 'LineWidth', 1.5);
% figure;plot(t(end-del:end),DDY(end-del:end),'r', 'LineWidth', 1.5);hold on;plot(t(end-del:end),D_qx1(end-del:end,2),'k.', 'LineWidth', 1.5);
% figure;plot(t(end-del:end),DDZ(end-del:end),'r', 'LineWidth', 1.5);hold on;plot(t(end-del:end),D_qx1(end-del:end,3),'k.', 'LineWidth', 1.5);
tic;
% Construct Library
[Theta,Sym]=LIBA(X,X_OrderMax,Trig_OrderMax,nonsmooth_OrderMax);
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
switch X_OrderMax
    case 3
        Xi_true=[-10	 4	0
            10	5	0
            0	0	-8/3
            0	0	0
            0	0	1
            0	-1	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0
            0	0	0];
    case 2
        Xi_true=[-10	4	0
            10	5	0
            0	0	-8/3
            0	0	0
            0	0	1
            0	-1	0
            0	0	0
            0	0	0
            0	0	0];
end

% switch X_OrderMax
%     case 3
%         Xi_true=[-10   -6	0
%             10	7.5	0
%             0	0	-8/3
%             0	0	0
%             0	0	1
%             0	-1	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0
%             0	0	0];
%     case 2
%         Xi_true=[-10 -6	0
%             10	7.5	0
%             0	0	-8/3
%             0	0	0
%             0	0	1
%             0	-1	0
%             0	0	0
%             0	0	0
%             0	0	0];
% end

RE = (Xi_optM(abs(Xi_true)>0) - Xi_true(abs(Xi_true)>0))./Xi_true(abs(Xi_true)>0);
MRE = max(abs(RE))*100;
fprintf('The maximum relative error %2.2f\n',MRE);

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
plot3(SM(1,:),SM(2,:),SM(3,:),'b--', 'LineWidth', 1.5);
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
scatter3(outputs.q, outputs.ObFu1, outputs.ObFu2, 30, distancesX, 'filled');
colorbar;
hold on;
scatter3(outputs.q(min_index_X), outputs.ObFu1(min_index_X), outputs.ObFu2(min_index_X), ...
    300, 'p', 'filled', 'MarkerEdgeColor',...
    'k', 'MarkerFaceColor', 'm');
% x_val = outputs.q(min_index_X);
% y_val = outputs.ObFu1(min_index_X);
% z_val = outputs.ObFu2(min_index_X);
% plot3([x_val, x_val], [y_val, y_val], [0, z_val], '--k', 'LineWidth', 2);
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

% Assuming you have already created a figure handle 'h'
% savefig(h, fullfile(savePath, ['FiguresFile_q', num2str(q1), '.fig']));
% 
% save(fullfile(savePath, ['data', num2str(q1), '.mat']));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%