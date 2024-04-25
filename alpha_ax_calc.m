% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function alpha = alpha_ax_calc(mq, rho_g,e_gm,sigma,rho_l)
% Function to determine the entrainment factor of axial symmetric flow (Milgram, 1983)

A = 7.598;
K = 0.165;
g = 9.81;

q = mq/ rho_g;
delta = e_gm; 
L_M = ( q^2 / ( g * delta^2 ) )^(1/5);
L_D = (sigma/ ( g *(rho_l-rho_g) ) )^0.5 / delta^(1/3);
F_B = delta^(2/5) * L_M / L_D;
alpha = K* F_B / ( A + F_B );

end


% LITERATURE:
% Milgram, J.H., 1983. Mean flow in round bubble plumes. Journal of Fluid Mechanics 133, 345. doi: 10.1017/S0 0221120830 01950