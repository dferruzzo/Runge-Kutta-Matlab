% ==================================================================
% Runge-Kutta Fehlberg (RK45) 
% autor: Diego Ferruzzo
% data: 26/10/2026
% ref: Numerical Methods for Engineers and Scientists Using MATLAB -
%      Estefandiari, Ramin, 2013
% ==================================================================
function [t,y,error] = rk45(f,t0,y0,h)
    k1 = f( t0, y0 );
    k2 = f( t0 + h/4, y0 + k1*h/4 );
    k3 = f( t0 + 3*h/8, y0 + 3*k1*h/32 + 9*k2*h/32 );
    k4 = f( t0 + 12*h/13, y0 + 1932*k1*h/2197 - 7200*k2*h/2197 + 7296*k3*h/2197 );
    k5 = f( t0 + h, y0 + 439*k1*h/216 - 8*k2*h + 3680*k3*h/513 - 845*k4*h/4104 );
    k6 = f( t0 + h/2, y0 - 8*k1*h/27 + 2*k2*h - 3544*k3*h/2565 + 1859*k4*h/4104 - 11*k5*h/40 );
    
    % RK4
    y = y0 + h*( 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5 );
    %y_rk5 = y0 + h*( 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55 );
    t = t0 + h;
    % estimativa do erro
    error = h*( k1/360 - 128*k3/4275 - 2197*k4/75240 + k5/50 + 2*k6/55 ); 
end