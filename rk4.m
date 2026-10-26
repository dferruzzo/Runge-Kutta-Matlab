% ==========================
% Runge Kutta 4a ordem
% autor: Diego Ferruzzo
% data: 26/10/2026
% ==========================
function [t, y] = rk4(f, t0, y0, h)
    k1 = f( t0, y0 );
    k2 = f( t0 + (h/2), y0 + (h*k1/2) );
    k3 = f( t0 + (h/2), y0 + (h*k2/2) );
    k4 = f( t0 + h, y0 + (h*k3) ); 
    y = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    t = t0 + h;
end