% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function [ d_p] = d_prim_calc( rho_g, rho_l, d_Orifice, sigma, q_Orifice,Temp )
% Function for calculating the primary bubble diameter according to
% Voit, Zeppenfeld, Mersmann 1987

%% Parameter
rho_g = rho_g;
rho_l = rho_l;
d_Orifice = d_Orifice;
sigma = sigma;
Temp = Temp;

%% Parameter calculation
drho = rho_l - rho_g;
eta_l = mu_Wasser_calc(Temp);

%% Constants
g = 9.81;

%% Initial conditions
d_Psigma= (6*d_Orifice * sigma / (drho*g ))^(1/3);        % Minimum bubble size

%% Iterative method for determining the primary bubble diameter
i = 1;
d_p(i) = d_Psigma;
dd_p = 1;

while abs(dd_p) > 0.00001
    i =  i+1;
    F_eta   = 15 * eta_l * q_Orifice / d_p(i-1);
    F_T     = 1.3 * rho_l * (q_Orifice / d_p(i-1))^2;
    F_sigma = pi * d_Orifice * sigma;
    d_p(i) = ( ( F_eta + F_T + F_sigma ) /(drho * g) * 6/pi )^(1/3);
    dd_p = (d_p(i)-d_p(i-1)) / d_p(i-1) ;
end

d_p = d_p(end);

end

% LITERATURE
% Voit, H., Zeppenfeld, R., Mersmann, A., 1987. Calculation of primary bubble vol- ume in gravitational and centrifugal fields. Chemical Engineering & Technology - CET 10, 99–103. doi: 10.1002/ceat.270100113 . http://doi.wiley.com/10.1002/ceat. 270100113 .
