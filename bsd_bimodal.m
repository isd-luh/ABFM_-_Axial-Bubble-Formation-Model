% ABFM - Axial Bubble Formation Model
% Copyright C 2024 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function [ z, a_Stutz, u_lzm, e_gm1, e_gm2, b,epsilon, v_mean1, v_mean2, f_a_Stutz ] = bsd_bimodal( a_Stutz, T, q_atm, p_atm, rho_l, sigma, Temp, d_Orifice, koal_cond, z_Auswertung )
% Function to obtain the bubble size distribution (BSD) from an axial flow
% The model bases on Bohne et al (2020).
% The balance equations for the mean bubble volume and gas fraction of the
% small and large bubble fraction base on Lehr et al. (2002)

%% Input Parameters with explainatation
T                   = T;                            % Water depth
q_atm               = q_atm;                        % Air volume flow at atmospheric condition
p_atm               = p_atm;                        % Atmospheric pressure;
Temp                = Temp;                         % Temperature    
rho_l               = rho_l;                        % Density of the fluid
sigma               = sigma;                        % Surface tension
d_Orifice           = d_Orifice;                    % Orifice Diameter
koal_cond           = koal_cond;                    % Coalescence Condition
z_Auswertung        = z_Auswertung;                 % Evaluation height for BSD    
a_Stutz             = a_Stutz;

% Constants for the axial symmetric flow
bc_coeff.alpha      = 0.08;     
bc_coeff.lambda     = 0.6;      
bc_coeff.u_rel1     = 0.3;    
bc_coeff.u_rel2     = 0.3;      
bc_coeff.ampli      = 1.0;

% Parameter fluid dynamics
lambda  = bc_coeff.lambda;                          % Ratio of eccentricity of liquid phase to gas phase      
u_rel1  = bc_coeff.u_rel1;                          % Bubble rise velocity Phase 1
u_rel2  = bc_coeff.u_rel2;                          % Bubble rise velocity Phase 2
ampli   = bc_coeff.ampli;                           % Amplification factor
alpha   = bc_coeff.alpha;                           % Entrainment factor

%% Constants
g = 9.81;                                           % gravity

%% Defining the coordinate system
%a_Stutz         = 0.0001:0.00001:0.5;
v_Stutz         = 4/3 * pi * a_Stutz.^3;

z_min           = d_Orifice * 6.2 ;
z_max           = T;        
z_Stutz         = linspace(z_min,z_max,1000);            

%% Determining the static pressure 
p_stat          = p_stat_calc(T-z_Stutz, rho_l, p_atm);                     % static pressure 

%% Determining the gas density 
rho_g           = rho_Luft_calc(Temp,p_stat);                               % density of air
drho_gdz        = ( rho_g(2) - rho_g(1)) / (z_Stutz(2)-z_Stutz(1));         % vertical gradient of density    

%% Determining the gas volume and mass flow at atmospheric condition
q_atm_Orifice   = q_atm ;
mq0             = q_atm_Orifice * rho_Luft_calc(Temp,p_atm);

%% Determinging the orifice flow 
[q_Orifice,M_Orifice,u_g_Orifice,d_prim,~,~] = orifice_flow_calc( q_atm_Orifice,d_Orifice,T,rho_l,p_atm,Temp,sigma );

%% Initial conditions
% Initial conditions of the fluid dynamics
M0      = M_Orifice + (2 * q_Orifice)  / u_g_Orifice * rho_l * g * z_min ;          % Momentum flow at starting point                          
e_gm0   = 0.5;                                                                      % Gas fraction at starting point    

e_gm10 = 0.01*e_gm0;
e_gm20 = 0.99*e_gm0;  

% Iterative procedure to determine the initial conditions of the axial flow b_0 and u_zm0
    k=1;
    u_lzm0(k)   = u_g_Orifice;
    b_0(k)      = d_Orifice/2;
    du_lzm0(k)  = 1;
    db_l0(k)    = 1;

while abs(db_l0(k)+du_lzm0(k)) > 0.0001
    k=k+1;
    
    b_0(k)      =  1/(lambda*((pi*rho_g(1)*(e_gm10*u_rel1 + e_gm20*u_rel2))/mq0 + (u_lzm0(k-1)*pi*rho_g(1)*(e_gm10 + e_gm20))/(mq0*(lambda^2 + 1)))^(1/2));
    u_lzm0(k)   = (2^(1/2)*M0^(1/2)*(2*lambda^2 + 1)^(1/2))/(ampli^(1/2)*b_0(k)*rho_l^(1/2)*pi^(1/2)*(2*lambda^2 - 2*e_gm20*lambda^2 - 2*e_gm10*lambda^2 + 1)^(1/2));

    du_lzm0(k)  = (u_lzm0(k)-u_lzm0(k-1))/u_lzm0(k-1);
    db_l0(k)    = ( b_0(k)-b_0(k-1) ) / b_0(k-1);
end

% Initial conditions of the bubble volumes
a_init     = d_prim/2*1.1447;           % Obtaining the initial bubble size; Rule: a_prim = Mode
v_mean20   = 4/3*pi*(a_init)^3 ;
v_mean10   = v_mean20/30; 

%% Euler forward integration procedure to solve the system of differential equations

% Defining step width Schrittweite
dz_max  = 0.00001;
z       = z_min : dz_max: z_max;
dz_vec  = diff(z);

% Defining Solution array 
y       = zeros(length(z),6);
y(1,:)  = [u_lzm0(end), e_gm10, e_gm20, b_0(end),v_mean10, v_mean20];

for i = 2:length(z)
    [dydz,epsilon(i),u_rel1_vec(i-1),u_rel2_vec(i-1)]   = bsd_bimodal_tgl( z(i-1), y(i-1,:), sigma, d_prim, rho_g, rho_l, drho_gdz, bc_coeff, koal_cond, z_Stutz);
    y(i,:)                                              = dydz'*dz_vec(i-1) + y(i-1,:) ;
end

u_lzm   = y(:,1);
e_gm1   = y(:,2);
e_gm2   = y(:,3);
b       = y(:,4);
v_mean1 = y(:,5);
v_mean2 = y(:,6);

%% Control segment
controlflag = false;
if controlflag == true
    u_rel1_vec(end+1) = u_rel1_vec(end);
    u_rel2_vec(end+1) = u_rel2_vec(end);

    mq = (b.^2.*e_gm1.*lambda^2*pi.*interp1(z_Stutz,rho_g,z').*(u_rel1_vec'*lambda^2 + u_lzm + u_rel1_vec'))./(lambda^2 + 1) + (b.^2.*e_gm2.*lambda^2*pi.*interp1(z_Stutz,rho_g,z').*(u_rel2_vec'*lambda^2 + u_lzm + u_rel2_vec'))./(lambda^2 + 1);
    figure; plot(mq,z)
    pause
end

%% Probability distribution
[~,i_Auswertung] = min( abs(z - z_Auswertung ) );      % Index for the position of the evaluation

v_tilde1    = v_Stutz./v_mean1(i_Auswertung);
f_v1        = e_gm1(i_Auswertung)/v_mean1(i_Auswertung)^2*(2/pi)^0.5 * 1./(3*v_tilde1) .* exp(-2/9.*log( exp(9/8).*v_tilde1 ).^2) ;

v_tilde2    = v_Stutz./v_mean2(i_Auswertung);
f_v2        = e_gm2(i_Auswertung)/v_mean2(i_Auswertung)^2*exp(-v_tilde2) ;

f_v         = f_v1 + f_v2;              % Resulting density distribution
f_v         = f_v/trapz(v_Stutz,f_v);   % PDF in volume-domain


% Transformation from volume to radius
f_a_Stutz    = 4 * pi * a_Stutz.^2 .* f_v;   % PDF in a-Raum                                  
% f_a_Stutz    = f_a_Stutz/trapz(a_Stutz,f_a_Stutz);

end
% Literature
% Tobias Bohne, Tanja Grieﬂmann , Raimund Rolfes, 2020. Development of and efficient buoyant jet 
% integral model of a bubble plume coupled with a population dynamics model for bubble breakup and 
% coalescence to predict the transmission loss of a bubble curtain. 
% International Journal of Multiphase Flow 132, 103436. doi: 10.1016/j.ijmultiphaseflow.2020.103436

% Lehr, F., Millies, M., Mewes, D., 2002. Bubble-Size distributions and flow fields in bubble columns.
%  AIChE Journal 48, 2426ñ2443. doi: 10.1002/aic.690481103 

