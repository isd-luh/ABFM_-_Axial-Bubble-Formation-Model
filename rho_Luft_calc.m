% # ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function rho_Luft = rho_Luft_calc(Temp,p)
% Function to determine the density of air

%% Material parameter of air from VDI Wärmeatlas
T_Luft_Stutz = [0, 25, 50];                                                          % Temperature in °C 
P_Luft_Stutz = [1*10^5, 5*10^5, 10*10^5];                                            % Pressure in bar
rho_Luft_Stutz = [1.2758, 1.1685, 1.0779;... 
                  6.3940, 5.8500, 5.3923;...
                  12.823, 11.717, 10.79       ];

%% Determine the density of air
rho_Luft =  interp2(T_Luft_Stutz,P_Luft_Stutz,rho_Luft_Stutz,Temp,p);

end

% LITERATURE:
% VDI Wärmeatlas 2006