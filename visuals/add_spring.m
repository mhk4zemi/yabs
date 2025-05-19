function add_spring(scale, offset, offset_axial)

    n_coils = 10;      % Number of coils
    amplitude = 1*scale;     % Amplitude of the spring
    length = 10;       % Total length of the spring
    points_per_coil = 100; % Points per coil for smoothness

    % Create the spring shape
    x = offset_axial + linspace(0, length, n_coils * points_per_coil);    % X-coordinates (length direction)
    y = offset + amplitude * sin(2 * pi * n_coils * x / length);   % Y-coordinates (sinusoidal shape)

    % Plot the spring
    figure(1);hold on
    plot(x, y, 'color',[0.6350 0.0780 0.1840], 'LineWidth', 2);  % 2D plot

end
