% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function [ u_t, du_tdv ] = bub_rise_calc( v,sigma,rho_l )
% Function to determine the rising velocity of a bubble in a non moving medium from Clift (1978) Eq. 7-3

g       = 9.81;
d_e     = (6*v/pi)^(1/3);

u_t     = (2.14*sigma / (rho_l * d_e) + 0.505 *g *d_e)^0.5;
du_tdv  = ((101*g)/(100*pi*((6*v)/pi)^(2/3)) - (107*sigma)/(25*rho_l*pi*((6*v)/pi)^(4/3)))/(2*((101*g*((6*v)/pi)^(1/3))/200 + (107*sigma)/(50*rho_l*((6*v)/pi)^(1/3)))^(1/2));

end

% Literature
% Clift, R. , Grace, J.R. , Weber, M.E. , 1978. Bubbles, Drops, and Particles. Academic Press, New York .