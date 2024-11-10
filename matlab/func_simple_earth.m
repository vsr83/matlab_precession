function k = func_simple_earth(t, y, mu)
% FUNC_LIB - Compute the time derivative of the degrees of freedom for 
% numerical integration including librations.
%
% INPUTS:
%   t          The time after epoch at evaluation.
%   y          The degrees of freedom at evaluation 
%              (6 x 1).
%   mu         The standard gravitational parameter for each point-mass
%              (au^3/d^2, num_targets)
%
% OUTPUTS:
%   k          The time derivative of the DoFs for numerical integration.
%
% REFERENCES:
%   [1]        Chen, Shen - New estimates of the inertia tensor and rotation 
%              of the triaxial nonrigid Earth, Journal of Geophysical
%              Research, 2010.
%   [2]        The IAU 2009 System of Astronomical Constants.
%   [3]        Folkner, Williams, Boggs, Park, Kuchynka - Planetary and 
%              Lunar Ephemerides DE430 and DE431, IPN Progress Report,
%              2014.

% Zonal harmonics J_2, J_3 and J_4 in Earth's potential [3].
Je = [0.00108262545, -0.00000253241, -0.000001616];

% Tesseral harmonics in Earth's potential.
CSnm = [];

% Extract Euler angles and their time derivatives from the DoFs.
phi    = y(1);%mod(y(1), 2 * pi);
phi1   = y(2);
theta  = y(3);%mod(y(3), 2 * pi);
theta1 = y(4);
psi    = y(5);%mod(y(5), 2 * pi);
psi1   = y(6);

T = 365.25;
%T = 1;
M = mod(t * 2 * pi / T, 2*pi);
r_point = matrix_rot1d(0) * [cos(M); sin(M); 0];
r_point = matrix_to_body(phi, theta, psi) * r_point;


%disp(r_point')

%disp(atan2d(r_point(2), r_point(1)))
%disp(asind(r_point(3) / norm(r_point)));
%disp([phi, theta, psi, r_point']);
%disp(r_point');

% Principal moments of inertia (kg m^2)  [1].
I_kgm2 = [8.0100829e37, 8.0102594e37, 8.0364807e37];

% Earth mass (kg) [2].
m_earth = 5.9722e24;

% Astronomical unit (m) [2].
au_meters = 149597870700;

% Principal moments of inertia per unit mass (au^2)
I = I_kgm2 / (m_earth * au_meters^2);

%I = [0.599310440525557e-9, 0.599323646168254e-9, 0.601285510864330e-9];

% Equatorial radius of Earth (au) [2].
a = 6378136.6 / au_meters;

mu_sun   = 2.959122082855911e-04;
[acc_point, T_body] = acc_body(r_point, a, mu_sun, Je, CSnm);

% Apply Euler equations to compute the angular accelerations of the Euler
% angles.
[phi2, theta2, psi2] = acc_euler_body(phi, theta, psi, phi1, theta1, psi1, T_body, I);

%disp([mod(t, 1), psi, psi2]);

k = zeros(6, 1);
% d/dt phi = phi1
k(1) = y(2);
% d/dt phi1 = phi2
k(2) = phi2;
% d/dt theta = theta1
k(3) = y(4);
% d/dt theta1 = theta2
k(4) = theta2;
% d/dt psi = psi1
k(5) = y(6);
% d/dt psi1 = psi2
k(6) = psi2;
