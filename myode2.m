% Load SynapseModel.m in pdetool and export the parameters
% Must export the variables b, p, e, t, c, a, f, d from pdetool
clf
% Determine length of simulation and time step
tlist = 0:0.01:0.5;
% Define initial conditions
u0 = zeros(size(p, 2), 1); % Set everything to zero
n = 11; % Number of nodes the initial concentration is spread over
Cs = 10; % Initial concentration
set = [181, 141, 140, 175, 105, 133, 176, 139, 104, 180, 138]; % Data points in vesicle
for i = 1:size(set, 2)
    u0(set(i)) = Cs / n;
end

% Define global sink term due to cleavage in the cleft for different conditions
a1 = '0';
a2 = '5';
c = '10';

% Define parameters specific to myasthenia gravis (MG)
mg_a = '5'; % Adjust this value based on the impact of MG on the sink term 'a'
mg_c = '25'; % Adjust this value based on the impact of MG on the sink term 'c'

% Solve PDE for normal condition and Sarin gas condition
u1 = parabolic(u0, tlist, b, p, e, t, c, a1, f, d);
u2 = parabolic(u0, tlist, b, p, e, t, c, a2, f, d);

% Solve PDE for MG condition
u_mg = parabolic(u0, tlist, b, p, e, t, mg_c, mg_a, f, d);

% Run multiple iterations for normal, Sarin gas, and MG conditions
for i = 1:3
    ui1 = u0 + u1(:, end);
    ui2 = u0 + u2(:, end);
    ui_mg = u0 + u_mg(:, end);
    
    uf1 = parabolic(ui1, tlist, b, p, e, t, c, a1, f, d);
    uf2 = parabolic(ui2, tlist, b, p, e, t, c, a2, f, d);
    uf_mg = parabolic(ui_mg, tlist, b, p, e, t, mg_c, mg_a, f, d);
    
    u1 = [u1, uf1];
    u2 = [u2, uf2];
    u_mg = [u_mg, uf_mg];
end

clear M
for j = 1:size(u1, 2)
    figure(1)
    subplot(1, 3, 1)
    pdesurf(p, t, u1(:, j))
    xlabel('x')
    ylabel('y')
    zlabel('ACh Concentration')
    title('Sarin Gas')
    colormap('hsv')
    caxis([0 .05])
    axis([1.5 4 0 10 0 .2])

    subplot(1, 3, 2)
    pdesurf(p, t, u2(:, j))
    xlabel('x')
    ylabel('y')
    zlabel('ACh Concentration')
    title('Normal')
    colormap('hsv')
    caxis([0 .05])
    axis([1.5 4 0 10 0 .2])

    subplot(1, 3, 3)
    pdesurf(p, t, u_mg(:, j))
    xlabel('x')
    ylabel('y')
    zlabel('ACh Concentration')
    title('Myasthenia Gravis')
    colormap('hsv')
    caxis([0 .05])
    axis([1.5 4 0 10 0 .2])

    % Calculate time in milliseconds for the current step
    % current_time_ms = tlist(j) * 1000;

    % Add a title with the current time step in milliseconds
    % sgtitle(['Time: ', num2str(current_time_ms), ' ms'])
    M(j) = getframe(figure(1));
end

movie2avi(M, 'moviecompressed', 'compression', 'Cinepak')
