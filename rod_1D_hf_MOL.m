clear variables
close all

% c(x,t,u,Du/Dx) * Du/Dt = x^(-m) * D(x^m * f(x,t,u,Du/Dx))/Dx + s(x,t,u,Du/Dx)
% 
%   m must be 0, 1, or 2, corresponding to slab, cylindrical, or spherical symmetry,
%   respectively. The coupling of the partial derivatives with respect to
%   time is restricted to multiplication by a diagonal matrix c(x,t,u,Du/Dx).
%   The diagonal elements of c are either identically zero or positive.
%   An entry that is identically zero corresponds to an elliptic equation and
%   otherwise to a parabolic equation. There must be at least one parabolic
%   equation. An entry of c corresponding to a parabolic equation is permitted
%   to vanish at isolated values of x provided they are included in the mesh
%   XMESH, and in particular, is always allowed to vanish at the ends of the
%   interval.

T_env = 25;

x = linspace(0,1,101);
t = linspace(0,1,250);

m = 0;
sol = pdepe(m,...
    @(x,t,u,dudx)heatrod(x,t,u,dudx,T_env),...
    @(x)heatic(x,T_env),...
    @(xl,ul,xr,ur,t)heatbc(xl,ul,xr,ur,t,T_env), ...
    x,t);

u = sol(:,:,1);

figure
plot(t,u(:,round(length(x)/2,0)));
xlabel('Time (s)')
ylabel('Temperature (°C)')
title('Temperature at the rod midpoint')

% figure
% xlabel('Time (s)')
% ylabel('Temperature (°C)')
% title('Rod temperature')
% my_ylims = [min(min(u))*0.9 max(max(u))*1.1];
% for i=1:length(t)
%     plot(x,u(i,:))
%     ylim(my_ylims); 
%     pause(0.1);
% end

%% Local functions
function [c,f,s] = heatrod(x,t,u,dudx,T_env)
c = 1; % Thermal capacity
f = dudx;

% T_env = 293; % Ambient temperature
h = 0.1; % Heat transfer coefficient (adjust as needed)
% Source term: Newton's law of cooling
s = -h * (u - T_env);
end
%----------------------------------------------
function u0 = heatic(x,T_env)
% n = 2.404825557695773;
u0 = T_env; % besselj(0,n*x);
end
%----------------------------------------------
function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t,T_env)
% Left boundary: Fixed temperature at 273 K
pl = ul - 100; 
ql = 0; 

% Right boundary: Fixed temperature at 273 K
pr = ur - 25;
qr = 0; 
end
%----------------------------------------------
