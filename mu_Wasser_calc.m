% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function mu = mu_Wasser_calc(T_Wasser)
% Function to determine the dynamic viscosity of Water 

%% Parameter
T_Wasser = T_Wasser;

%% Definition of supporting points
mu_Stutz        = [1792*10^(-6);890*10^(-6)];                                       % Dynamic viscosity of water from VDI-Wärmeatlas 2006                                             
T_Wasser_Stutz  = [0 ; 25];                                                         % Support points water temperature

%% Calculating the dynamic viscosity
mu = interp1(T_Wasser_Stutz,mu_Stutz,T_Wasser);

end

% Literatur
% VDI-Wärmeatlas 2006