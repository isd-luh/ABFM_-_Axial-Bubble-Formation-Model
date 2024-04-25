% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 


function [ dydz,epsilon,u_rel1,u_rel2 ] = bsd_bimodal_tgl( z, y, sigma, d_prim, rho_g_Stutz, rho_l, drho_gdz, bc_coeff, koal_cond, z_Stutz   )
% Function to calculate the derivative for each step size of the euler forward integration procedure

%% Field values
u_lzm   = y(1);
e_gm1   = y(2);
e_gm2   = y(3);
b       = y(4);
v_mean1 = y(5);
v_mean2 = y(6);

%% constants
g           = 9.81;                % gravity

%% Phase control
if e_gm2/(e_gm1+e_gm2) < 0.0001 || v_mean2 < v_mean1
    ph1act = 1;
    ph2act = 0;
else
    ph1act = 1;
    ph2act = 1;
end
%% Interpolating the density on integration points
rho_g       = interp1(z_Stutz,rho_g_Stutz,z,'linear');

%% Parameter
lambda  = bc_coeff.lambda;    % Relation between the radial excentricity of fluid phase and gas phase      
ampli   = bc_coeff.ampli;     % Amplification factor

%% Algebraic functions

% bubble rising velocity
[ u_rel1, du_rel1dvmean1 ] = bub_rise_calc( v_mean1,sigma,rho_l );
[ u_rel2, du_rel2dvmean2 ] = bub_rise_calc( v_mean2,sigma,rho_l );
% u_rel1=0.3;
% u_rel2=0.3;
% du_rel1dvmean1=0;
% du_rel2dvmean2=0;

if z < d_prim/2
    epsilon     = 0;
    u_krit      = 0;
else
    epsilon     = -(g*lambda^2*(e_gm1*u_rel1 + e_gm2*u_rel2)*(exp(-2/lambda^2) - 1))/2 *koal_cond(2) ;
    %epsilon 	= Diss_calc(u_lzm, b)*koal_cond(2);
    u_krit      = 0.08*koal_cond(1);
end

u_res           = u_krit;
r_interm        = b * 2^0.5;          % Intermittenz zone

%% Obtaining the fluid dynamic parameters
mq      = pi*b^2*rho_g*(e_gm1*u_lzm + e_gm2*u_lzm + e_gm1*lambda^2*u_rel1 + e_gm2*lambda^2*u_rel2) - (pi*b^2*u_lzm*rho_g*(e_gm1 + e_gm2))/(lambda^2 + 1);  % Massenstrom
alpha   = alpha_ax_calc(mq, rho_g,(e_gm1+e_gm2) ,sigma, rho_l);

%% Coefficients
A(1,1) = -(ampli*b^2*rho_l*u_lzm*pi*(2*e_gm1*lambda^2 + 2*e_gm2*lambda^2 - 2*lambda^2 - 1))/(2*lambda^2 + 1);
A(1,2) = -(pi*ampli*b^2*lambda^2*rho_l*u_lzm^2)/(2*lambda^2 + 1);
A(1,3) = -(pi*ampli*b^2*lambda^2*rho_l*u_lzm^2)/(2*lambda^2 + 1);
A(1,4) = -(ampli*b*rho_l*u_lzm^2*pi*(2*e_gm1*lambda^2 + 2*e_gm2*lambda^2 - 2*lambda^2 - 1))/(2*lambda^2 + 1);
A(1,5) = 0;
A(1,6) = 0;
 
A(2,1) = -(b^2*rho_l*pi*(e_gm1*lambda^2 + e_gm2*lambda^2 - lambda^2 - 1))/(lambda^2 + 1);
A(2,2) = -(pi*b^2*lambda^2*rho_l*u_lzm)/(lambda^2 + 1);
A(2,3) =-(pi*b^2*lambda^2*rho_l*u_lzm)/(lambda^2 + 1);
A(2,4) = -(2*b*rho_l*u_lzm*pi*(e_gm1*lambda^2 + e_gm2*lambda^2 - lambda^2 - 1))/(lambda^2 + 1);
A(2,5) = 0;
A(2,6) = 0;
 
A(3,1) = (pi*b^2*e_gm1*lambda^2*rho_g)/(lambda^2 + 1);
A(3,2) = pi*b^2*lambda^2*u_rel1*rho_g + (pi*b^2*lambda^2*u_lzm*rho_g)/(lambda^2 + 1);
A(3,3) = 0;
A(3,4) = 2*pi*b*e_gm1*lambda^2*u_rel1*rho_g + (2*pi*b*e_gm1*lambda^2*u_lzm*rho_g)/(lambda^2 + 1);
A(3,5) = pi*b^2*e_gm1*lambda^2*rho_g*du_rel1dvmean1;
A(3,6) = 0;
 
A(4,1) = (pi*b^2*e_gm2*lambda^2*rho_g)/(lambda^2 + 1);
A(4,2) =  0;
A(4,3) = pi*b^2*lambda^2*u_rel2*rho_g + (pi*b^2*lambda^2*u_lzm*rho_g)/(lambda^2 + 1);
A(4,4) = 2*pi*b*e_gm2*lambda^2*u_rel2*rho_g + (2*pi*b*e_gm2*lambda^2*u_lzm*rho_g)/(lambda^2 + 1);
A(4,5) = 0;
A(4,6) = pi*b^2*e_gm2*lambda^2*rho_g*du_rel2dvmean2;

A(5,1) = (pi*b^2*e_gm1*lambda^2*v_mean1*rho_g)/(lambda^2 + 1);
A(5,2) = pi*b^2*lambda^2*v_mean1*u_rel1*rho_g + (pi*b^2*lambda^2*u_lzm*v_mean1*rho_g)/(lambda^2 + 1);
A(5,3) = 0;
A(5,4) = 2*pi*b*e_gm1*lambda^2*v_mean1*u_rel1*rho_g + (2*pi*b*e_gm1*lambda^2*u_lzm*v_mean1*rho_g)/(lambda^2 + 1);
A(5,5) = pi*b^2*e_gm1*lambda^2*u_rel1*rho_g + (pi*b^2*e_gm1*lambda^2*u_lzm*rho_g)/(lambda^2 + 1) + pi*b^2*e_gm1*lambda^2*v_mean1*rho_g*du_rel1dvmean1;
A(5,6) = 0;

A(6,1) = (pi*b^2*e_gm2*lambda^2*v_mean2*rho_g)/(lambda^2 + 1);
A(6,2) = 0;
A(6,3) = pi*b^2*lambda^2*v_mean2*u_rel2*rho_g + (pi*b^2*lambda^2*u_lzm*v_mean2*rho_g)/(lambda^2 + 1);
A(6,4) = 2*pi*b*e_gm2*lambda^2*v_mean2*u_rel2*rho_g + (2*pi*b*e_gm2*lambda^2*u_lzm*v_mean2*rho_g)/(lambda^2 + 1);
A(6,5) = 0;
A(6,6) = pi*b^2*e_gm2*lambda^2*u_rel2*rho_g + (pi*b^2*e_gm2*lambda^2*u_lzm*rho_g)/(lambda^2 + 1) + pi*b^2*e_gm2*lambda^2*v_mean2*rho_g*du_rel2dvmean2;
 
%% R.H.S
f(1) = pi*b^2*g*lambda^2*rho_l*(e_gm1 + e_gm2);
f(2) = 2*pi*alpha*b*rho_l*u_lzm;
f(3) = (31043*b^2*e_gm1^2*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1)*((2150310427208497*v_mean1)/1125899906842624)^(2/3))/(20000*v_mean2) - (pi*b^2*e_gm1*lambda^2*u_lzm*drho_gdz)/(lambda^2 + 1) - (509*b^2*e_gm2*lambda^2*ph1act*ph2act*pi*rho_g*epsilon^(1/15)*(sigma/rho_l)^(2/5)*(exp(-r_interm^2/(b^2*lambda^2)) - 1))/(2000*v_mean2^(4/9)) - (b^2*e_gm1*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*exp(-(2*r_interm^2)/(b^2*lambda^2))*rho_g*(31043*e_gm1*((2150310427208497*v_mean1)/1125899906842624)^(2/3) + 2256*e_gm2*((2150310427208497*v_mean2)/1125899906842624)^(2/3) + 2256*e_gm2*((10751552136042485*v_mean1)/1125899906842624)^(2/3) + 4512*e_gm2*((10751552136042485*v_mean1)/1125899906842624)^(1/3)*((2150310427208497*v_mean2)/1125899906842624)^(1/3)))/(20000*v_mean2) - pi*b^2*e_gm1*lambda^2*u_rel1*drho_gdz + (141*b^2*e_gm1*e_gm2*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(((10751552136042485*v_mean1)/1125899906842624)^(1/3) + ((2150310427208497*v_mean2)/1125899906842624)^(1/3))^2*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1))/(1250*v_mean2);
f(4) = (509*b^2*e_gm2*lambda^2*ph1act*ph2act*pi*rho_g*epsilon^(1/15)*(sigma/rho_l)^(2/5)*(exp(-r_interm^2/(b^2*lambda^2)) - 1))/(2000*v_mean2^(4/9)) - (pi*b^2*e_gm2*lambda^2*u_lzm*drho_gdz)/(lambda^2 + 1) - pi*b^2*e_gm2*lambda^2*u_rel2*drho_gdz + (b^2*e_gm1*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*exp(-(2*r_interm^2)/(b^2*lambda^2))*rho_g*(31043*e_gm1*((2150310427208497*v_mean1)/1125899906842624)^(2/3) + 2256*e_gm2*((2150310427208497*v_mean2)/1125899906842624)^(2/3) + 2256*e_gm2*((10751552136042485*v_mean1)/1125899906842624)^(2/3) + 4512*e_gm2*((10751552136042485*v_mean1)/1125899906842624)^(1/3)*((2150310427208497*v_mean2)/1125899906842624)^(1/3)))/(20000*v_mean2) - (31043*b^2*e_gm1^2*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1)*((2150310427208497*v_mean1)/1125899906842624)^(2/3))/(20000*v_mean2) - (141*b^2*e_gm1*e_gm2*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(((10751552136042485*v_mean1)/1125899906842624)^(1/3) + ((2150310427208497*v_mean2)/1125899906842624)^(1/3))^2*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1))/(1250*v_mean2);
f(5) = (3041*b^2*e_gm1*lambda^2*ph1act*v_mean1^(14/9)*pi*exp(-(1713*(sigma/rho_l)^(9/10))/(2000*v_mean1^(1/2)*epsilon^(3/5)))*rho_g*epsilon^(19/15)*(rho_l/sigma)^(7/5)*(exp(-r_interm^2/(b^2*lambda^2)) - 1))/5000 - (pi*b^2*e_gm1*lambda^2*u_lzm*v_mean1*drho_gdz)/(lambda^2 + 1) - (3463*b^2*e_gm1^2*lambda^2*ph1act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1)*((2150310427208497*v_mean1)/1125899906842624)^(2/3))/20000 - (b^2*e_gm1*lambda^2*ph1act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*exp(-(2*r_interm^2)/(b^2*lambda^2))*rho_g*(31043*e_gm1*ph2act*v_mean1*((2150310427208497*v_mean1)/1125899906842624)^(2/3) - 3463*e_gm1*v_mean2*((2150310427208497*v_mean1)/1125899906842624)^(2/3) + 2256*e_gm2*ph2act*v_mean1*((2150310427208497*v_mean2)/1125899906842624)^(2/3) + 2256*e_gm2*ph2act*v_mean1*((10751552136042485*v_mean1)/1125899906842624)^(2/3) + 4512*e_gm2*ph2act*v_mean1*((10751552136042485*v_mean1)/1125899906842624)^(1/3)*((2150310427208497*v_mean2)/1125899906842624)^(1/3)))/(20000*v_mean2) - pi*b^2*e_gm1*lambda^2*v_mean1*u_rel1*drho_gdz - (509*b^2*e_gm2*lambda^2*ph1act*ph2act*v_mean1*pi*rho_g*epsilon^(1/15)*(sigma/rho_l)^(2/5)*(exp(-r_interm^2/(b^2*lambda^2)) - 1))/(2000*v_mean2^(4/9)) + (31043*b^2*e_gm1^2*lambda^2*ph1act*ph2act*u_res*v_mean1*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1)*((2150310427208497*v_mean1)/1125899906842624)^(2/3))/(20000*v_mean2) + (141*b^2*e_gm1*e_gm2*lambda^2*ph1act*ph2act*u_res*v_mean1*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(((10751552136042485*v_mean1)/1125899906842624)^(1/3) + ((2150310427208497*v_mean2)/1125899906842624)^(1/3))^2*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1))/(1250*v_mean2);
f(6) = (b^2*lambda^2*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*exp(-(2*r_interm^2)/(b^2*lambda^2))*rho_g*(4250*e_gm2^2*((2150310427208497*v_mean2)/1125899906842624)^(2/3) + 2256*e_gm1*e_gm2*((10751552136042485*v_mean1)/1125899906842624)^(2/3) + 2256*e_gm1*e_gm2*((2150310427208497*v_mean2)/1125899906842624)^(2/3) + 31043*e_gm1^2*ph1act*((2150310427208497*v_mean1)/1125899906842624)^(2/3) + 4512*e_gm1*e_gm2*((10751552136042485*v_mean1)/1125899906842624)^(1/3)*((2150310427208497*v_mean2)/1125899906842624)^(1/3) + 2256*e_gm1*e_gm2*ph1act*((10751552136042485*v_mean1)/1125899906842624)^(2/3) + 2256*e_gm1*e_gm2*ph1act*((2150310427208497*v_mean2)/1125899906842624)^(2/3) + 4512*e_gm1*e_gm2*ph1act*((10751552136042485*v_mean1)/1125899906842624)^(1/3)*((2150310427208497*v_mean2)/1125899906842624)^(1/3)))/20000 - (pi*b^2*e_gm2*lambda^2*u_lzm*v_mean2*drho_gdz)/(lambda^2 + 1) - pi*b^2*e_gm2*lambda^2*v_mean2*u_rel2*drho_gdz + (509*b^2*e_gm2*lambda^2*ph2act*v_mean2^(5/9)*pi*rho_g*epsilon^(1/15)*(sigma/rho_l)^(2/5)*(exp(-r_interm^2/(b^2*lambda^2)) - 1))/2000 - (17*b^2*e_gm2^2*lambda^2*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1)*((2150310427208497*v_mean2)/1125899906842624)^(2/3))/80 - (141*b^2*e_gm1*e_gm2*lambda^2*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(((10751552136042485*v_mean1)/1125899906842624)^(1/3) + ((2150310427208497*v_mean2)/1125899906842624)^(1/3))^2*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1))/1250 + (509*b^2*e_gm2*lambda^2*ph1act*ph2act*v_mean2^(5/9)*pi*rho_g*epsilon^(1/15)*(sigma/rho_l)^(2/5)*(exp(-r_interm^2/(b^2*lambda^2)) - 1))/2000 - (31043*b^2*e_gm1^2*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1)*((2150310427208497*v_mean1)/1125899906842624)^(2/3))/20000 - (141*b^2*e_gm1*e_gm2*lambda^2*ph1act*ph2act*u_res*pi^2*exp(-(1899241518582531/(2251799813685248*(e_gm1 + e_gm2)^(1/3)) - 1)^2)*rho_g*(((10751552136042485*v_mean1)/1125899906842624)^(1/3) + ((2150310427208497*v_mean2)/1125899906842624)^(1/3))^2*(exp(-(2*r_interm^2)/(b^2*lambda^2)) - 1))/1250;

%% System of equation
dydz = inv(A) * f';


end

