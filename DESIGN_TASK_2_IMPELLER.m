clear all

% Number of points
n = 100

% Initialization of variables
r = linspace(0.1, 0.5, n)
rho = zeros(1, n)
b_n = zeros(1, n)
h = zeros(1, n)

% Initial data

PR = 3                                  % Pressure ratio
m_flow = 2                              % Mass flow
omega = 1000 * (2*pi)/60                % Angular velocity 
z = 50                                  % Number of blades
b_n(1) = 0.1                            % Initial blade height
eta = 0.75                              % Efficiency

k = 1.4                                 % Ratio of specific heats (gamma)
T(1) = 309.34                           % Initial temperature
rho_initial = 1.28
rho(1) = rho_initial                    % Initial air density
R = 287                                 % Perfect gas constant for air (in J/(kg.K))
P_0 = 1*10^5                            % Initial stagnation pressure

r = linspace(0.1, 0.5, 100)

sigma = 1 - (0.63*pi)/z                 % Slip factor
c_p = (k*R)/(k-1)                       % Specific hear for cold air (kJ/kg)

% Calculation of radial velocity (constant throughout the compressor)
C_r = m_flow / (2*pi*r(1)*b_n(1)*rho(1))

% Calculation of tagential velocity
C_w = sigma * omega .* r

% Calculation of total velocity
C = sqrt(C_w.^2 + C_r^2)

% Entalpy at the inlet and outlet of the impeller
T(n) = T(1)*PR^((k-1)/k)
h(1) = c_p*T(1)
h(n) = c_p*T(n) % Valid for isentropic evolution

h(n) = h(1) + (h(n)-h(1))/eta

% Entalpy as a function of r
m = (h(n)-h(1))/(r(n)-r(1))
b = h(1) - m*r(1)

h = m.*r + b

% Since h = c_p * T => T(r) = h/c_p
T = h./c_p

% With static temperature and total velocity, we can get stagnation
% temperature, which will be used for the pressure calculation
T_0 = T + (C.^2)/(2*c_p) 

% Static pressure from stagnation values
P = P_0 * (T./T_0).^(k/(k-1))

% Density along the radius
rho = P./(R.*T)
rho = rho + (rho_initial-rho(1))

% Blade height
b_n = m_flow ./ (2*pi.*r*C_r.*rho)

% Plots

subplot(3,3,1);
plot(r, b_n); xlabel('Radius r (m)'); ylabel('Blade Height b_n (m)'); 
grid on;

subplot(3,3,2);
plot(r, h./(10^3)); xlabel('Radius r (m)'); ylabel('Enthalpy h (kJ/kg)');
grid on;

subplot(3,3,3);
plot(r, P./(10^5)); xlabel('Radius r (m)'); ylabel('Static Pressure P (bar)'); 
grid on;

subplot(3,3,4);
plot(r, rho); xlabel('Radius r (m)'); ylabel('Density rho (kg/m^3)'); 
grid on;

subplot(3,3,5);
plot(r, T); xlabel('Radius r (m)'); ylabel('Static Temperature T (K)'); 
grid on;

subplot(3,3,6);
plot(r, C); xlabel('Radius r (m)'); ylabel('Total Velocity C (m/s)'); 
grid on;

subplot(3,3,7);
plot(r, C_w); xlabel('Radius r (m)'); ylabel('Tangential Velocity C_w (m/s)'); 
grid on;

sgtitle('Impeller Flow Properties');