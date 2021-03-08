function [Bold] = BOLD_MFM(NA)

% ----------------------------------------------------------------------- %
% <<<<<<<<<<<<<<< complementary code for MFM algorithm >>>>>>>>>>>>>>>>>> %
% ----------------------------------------------------------------------- %

global itaus itauf itauo ialpha Eo dt
dt  = 0.01;
n_t=length(NA);
t_min = 10;
n_min = round(t_min/dt)+1;
taus   = 0.65;
tauf   = 0.41;
tauo   = 0.98;
alpha  = 0.32;
itaus  = 1/taus;
itauf  = 1/tauf;
itauo  = 1/tauo;
ialpha = 1/alpha;
Eo     = 0.34;
vo     = 0.02;
k1     = 7*Eo;
k2     = 2;
k3     = 2*Eo-0.2;
x0  = [0 1 1 1];
x = zeros(n_t,4);
x(1,:) = x0;
for n = 1:n_t-1
    x(n+1,1) = x(n,1) + dt*( NA(n)-itaus*x(n,1)-itauf*(x(n,2)-1) );
    x(n+1,2) = x(n,2) + dt*x(n,1);
    x(n+1,3) = x(n,3) + dt*itauo*(x(n,2)-x(n,3)^ialpha);
    x(n+1,4) = x(n,4) + dt*itauo*(x(n,2)*(1-(1-Eo)^(1/x(n,2)))/Eo - (x(n,3)^ialpha)*x(n,4)/x(n,3));
end
v  = x(n_min:end,3);
q  = x(n_min:end,4);
Bold  = 100/Eo*vo*( k1.*(1-q) + k2*(1-q./v) + k3*(1-v) );
clear x;