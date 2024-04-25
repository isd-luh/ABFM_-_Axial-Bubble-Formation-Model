% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function d_Brauer = d_Brauer_calc(Q_orifice, D_orifice)
% Function for determining the secondary bubble size according Brauer as given in
% VDI Wärmeatlas (2006)

g = 9.81 ;
A_orifice = pi * D_orifice^2 / 4 ;
u_orifice = Q_orifice / A_orifice ;

%% Determining Froude number
Fr_N = u_orifice^2 / ( D_orifice * g ) ;                                                    % Froude Zahl der Partikelbildung

%% Determining the bubble diameter
d_Brauer = D_orifice * 0.72 * Fr_N^(1/6) ;                                                  % Secondary Bubble size after jet decay                

end

% Literature
% VDI-Wärmeatlas 2006