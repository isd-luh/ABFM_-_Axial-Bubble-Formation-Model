% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function p_stat = p_stat_calc(T, rho_w, p_atm)
% Function to determine the static pressure in water 

%% Parameter
g = 9.81;                                   % gravity

%% determine the static pressure in water
p_stat = rho_w* g * abs(T)+ p_atm;

end