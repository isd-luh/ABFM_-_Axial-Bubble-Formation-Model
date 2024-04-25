% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function [ q_Orifice,M_Orifice,u_g_Orifice,d_prim,d_sek,We_N ] = orifice_flow_calc( q_atm,d_Orifice,T,rho_l,p_atm,Temp,sigma )
% Function to determine the air flow, out flow velocity, momentum flow,
% secondary bubble size and primary bubble size of a nozzle flow

%% Parameter
q_atm = q_atm;
d_Orifice= d_Orifice;
T = T;
rho_l = rho_l;
p_atm = p_atm;
Temp = Temp;
sigma = sigma;

%% Constants
g = 9.81;

%% Berechnung Düsenströmung
mq              = rho_Luft_calc(Temp, p_atm) * q_atm;                                % air mass flow
p_stat_Orifice  = p_stat_calc(T, rho_l, p_atm);                                      % static pressure
rho_g_Orifice   = rho_Luft_calc(Temp, p_stat_Orifice);                               % air density
A_Orifice       = pi*(d_Orifice/2)^2;                                                % cross section nozzle
u_g_Orifice     = mq / ( A_Orifice * rho_g_Orifice );                                % flow veloctiy in nozzle
M_Orifice       = mq * u_g_Orifice;                                                  % momentum flow through nozzle                 
q_Orifice       = mq / rho_g_Orifice;                                                % air volume flow   

%% Blasengröße
d_prim  = d_prim_calc( rho_g_Orifice, rho_l, d_Orifice, sigma, q_Orifice,Temp  );    % Primary bubble size 
d_sek   = d_Brauer_calc(q_Orifice, d_Orifice);                                       % Secondary bubble size 

%% Dimensionslose Kennzahlen
% VDI Wärmeatlas 
We_N = u_g_Orifice^2 * d_Orifice * rho_g_Orifice / sigma;                           % Weber number of particle size

end

% LITERATURE:
% VDI Wärmeatlas 