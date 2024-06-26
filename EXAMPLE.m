% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% EXAMPLE - Script to start the calculation of the bubble size distribution in an axial flow above a nozzle 
% The Axial Bubble Formation Model bases on Bohne et al (2020).
% Red curve in Fig. 11 in Bohne et al 2020

%% Input-Parameter 
T                   = 1;                            % Water depth (m)
q_atm_Orifice       = 2.2516*10^-4;                        % Air volume flow at atmospheric pressure through nozzle (m^3/s)
p_atm               = 101325;                       % Atmospheric pressure (Pa);
Temp                = 20;                           % Temperature (�C)
rho_l               = 1000;                         % Water density (kg/m^3)
sigma               = 0.072;                        % Surface tension (N/m)
d_Orifice           = 0.002;                        % Nozzle diameter (m)      
koal_cond           = [1,1];                        % Coalescence conditions: [1,1]: Break up and coalescence, [1,0]: Coalescence and no Break up, [0,1]: Break up and no Coalescence 
z_Auswertung        = 0.55;                            % Evaluation height above nozzle (m)
a_Klassen_Mitte     = 0.0001:0.0001:0.05;           % Bubble radius supporting points (m)

%% Solution of the fluid dynamics and transport equations for an axial flow 

[ z, a_Stutz, ~, ~, ~, ~,~, ~, ~, f_a_Stutz ] = bsd_bimodal( a_Klassen_Mitte, T, q_atm_Orifice, p_atm, rho_l, sigma, Temp, d_Orifice, koal_cond, z_Auswertung );

%% PDF of the gas volume
da          = a_Stutz(2)-a_Stutz(1);
V_i         = f_a_Stutz .* 4/3 * pi .* a_Stutz.^3 ;
V           = sum(V_i*da);
h_V         = V_i/V ;

%% Output
% figure;plot(a_Stutz,f_a_Stutz)
figure;plot(a_Stutz,h_V)


% Literature
% Tobias Bohne, Tanja Grie�mann , Raimund Rolfes, 2020. Development of and efficient buoyant jet integral model of a bubble plume coupled with a population dynamics model for bubble breakup and coalescence to predict the transmission loss of a bubble curtain. International Journal of Multiphase Flow 132, 103436. doi: 10.1016/j.ijmultiphaseflow.2020.103436


