% ======================================================
% Arquivo: Test - Runge-Kutta 4
% Autor: Diego Ferruzzo
% Data: 26/10/2026
%
%   Simula o atrator de Lorenz utilizando a funçao rk4.m
%
% ======================================================
%
clear;
close all;
clc;
%% Definiçao da funçao de Lorenz
function dx = lorenz(x, params)
    sigma = params(1);
    rho = params(2);
    beta = params(3);
    dx(1) = sigma*( x(2) - x(1) );
    dx(2) = x(1)*( rho - x(3) ) - x(2);
    dx(3) = x(1)*x(2) - beta*x(3);
end

params(1) = 10; % sigma
params(2) = 28; % rho
params(3) = 8/3; % beta
odefun = @(t,x) lorenz(x, params);
%% Integraçao numerica utilizando RK4 e RK45
t0 = 0; % tempo inicial
tf = 100; % tempo final
y0 = [0.9,0,0]; % condiçao inicial
dt = 1e-3; % passo de integraçao

N = floor((tf-t0)/dt); % numero de passos
t_rk4 = nan(N,1); % vetor de tempo vazio
y_rk4 = nan(N,3); % vetor de saida vazio

t_rk45 = nan(N,1); % vetor de tempo vazio
y_rk45 = nan(N,3); % vetor de saida vazio
error = zeros(N,3); % vetor de erro estimado

% Condiçao inicial no primeiro elemento da saida
t_rk4(1) = t0;
t_rk45(1) = t0;

y_rk4(1,:) = y0;
y_rk45(1,:) = y0;

% integraçao numerica
for i = 1:N-1
    [t_rk4(i+1), y_rk4(i+1,:)] = rk4(odefun, t_rk4(i), y_rk4(i,:), dt);
    [t_rk45(i+1), y_rk45(i+1,:), error(i+1,:)] = rk45(odefun, t_rk45(i), y_rk45(i,:), dt);
end
%% Graficos
figure(1);
subplot(1,2,1);
plot(y_rk4(:,1),y_rk4(:,2),'b-');
hold on;
plot(y_rk45(:,1),y_rk45(:,2),'r--');
xlabel('x');
ylabel('y');
grid;
subplot(1,2,2);
plot(y_rk4(:,2),y_rk4(:,3),'b-');
hold on;
plot(y_rk45(:,2),y_rk45(:,3),'r--');
xlabel('x');
ylabel('z');
grid;
%% processando o erro
% Calculate mean and standard deviation across experiments
error_norm = vecnorm(error,2,2);
media = mean(error, 2);
desvio = std(error, 0, 2);

% Plot mean with standard deviation bars
figure(2);
errorbar(t_rk45, media, desvio);
grid on;
xlabel('t');
ylabel('Mean \pm standard deviation');
legend('Mean \pm std');
%%

figure(3);
hold on;

fill([t_rk45; flipud(t_rk45)], ...
     [media + desvio; flipud(media - desvio)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.5);

plot(t_rk45, media, 'k', 'LineWidth', 1.5);

grid on;
xlabel('t');
ylabel('Mean \pm standard deviation');
legend('Standard deviation', 'Mean');
hold off;
