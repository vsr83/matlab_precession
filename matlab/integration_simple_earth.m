phi = 0.0001;
phi1 = 0;
theta = deg2rad(23.4);
theta1 = 0;
psi = 0.0001;
psi1 = 2 * pi;

y_initial = [phi, phi1, theta, theta1, psi, psi1];

% Time step size in days.
JT_timestep = 0.02;
% Start date 1969-Jun-28 00:00:00.0000
JT_epoch = 2440400.50;
% Compute number of timesteps.
num_timesteps = floor(365.25 * 10 / JT_timestep + 1);

Y = zeros(num_timesteps, 6);

y = y_initial';
t = 0;

for timestep = 1:num_timesteps
    if mod(timestep, 10000) == 0
        disp([timestep, num_timesteps])
    end

    JT = t + JT_epoch;
    [y, t] = runge4(@func_simple_earth, t, y, JT_timestep, 1);

    Y(timestep, :) = y';
end