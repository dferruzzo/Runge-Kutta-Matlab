% =============================================================
% Arquivo: Teste da funçao rk45_h_adj.m
% Autor: Diego Ferruzzo
% Data: 26/10/2026
%
%   Simula o atrator de Lorenz utilizando a funçao rk45_h_adj.m
%
% =============================================================
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
%% 
t0 = 0; % tempo inicial
tf = 100; % tempo final
y0 = [0.9,0,0]; % condiçao inicial
h_min = 1e-4; % passo de integraçao minimo
h_max = 1e-2; % passo de integraçao maximo
tol = 1e-9;   % tolerancia
h = h_max;
tempo = t0;

N = floor((tf-t0)/h_min);
t_hist = NaN(N,1);
y_hist = NaN(N,3);
error_hist = NaN(N,3);
h_hist = NaN(N,1);

i=1;

y_hist(1,:) = y0;
t_hist(1,:) = t0;
error_hist(1,:) = [0,0,0];
h_hist(1,1) = h;

t = t0;

while t < tf
    [t, y, error, h_adj] = rk45_h_adj(odefun, t0, y0, h, tol);
    if h_adj < h_min
        h_adj = h_min;
    end
    if h_adj > h_max
        h_adj = h_max;
    end
    
    t0 = t;
    y0 = y;
    h = h_adj;

    y_hist(i+1,:) = y;    
    t_hist(i+1,:) = t;
    error_hist(i+1,:) = error;
    h_hist(i+1,:) = h;

    i = i + 1;
end

display('Done');

%% Remove all NaN
M  = i-1;
y_clean = y_hist(1:M,:);
t_clean = t_hist(1:M);
h_clean =  h_hist(1:M);
error_clean = error_hist(1:M,:); 


%% plots
figure(1);
plot(t_clean,h_clean);
hold on;
yline(h_min, 'b');
yline(h_max, 'r');
title('Passo de integração');
xlabel('Tempo (s)');
ylabel('Passo (s)');

figure(2);
plot(y_clean(:,1),y_clean(:,2));
title('Atrator de Lorenz');
xlabel('x');
ylabel('y');

figure(3);
plot(y_clean(:,2),y_clean(:,3));
title('Atrator de Lorenz');
xlabel('y');
ylabel('z');

figure(3);
plot(t_clean, error_clean);
title('Erro estimado');
xlabel('Tempo (s)');
ylabel('Erro estimado');
