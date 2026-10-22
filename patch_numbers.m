freq = 4.9e9; % resonant frequency
c = 3e8; % speed of light in a vacuum
eps_r = 4.3; % relative permittivity
R_in_y0 = 50; % desired input resistance
Z_0 = 50; % desired input impedance

W = (c / (2*freq)) * sqrt(2 / (eps_r+1)) % width of patch

h_values = [0.8, 1.6, 2.4] * 1e-3; % candidate substrate heights

for i = 1:length(h_values)
    h = h_values(i); % substrate height

    eps_reff = ((eps_r+1) / 2) + ((eps_r-1)/2)*(1 + 12*(h/W))^(-1/2); % effective dialectric constant

    delta_L = 0.412 * h * ((eps_reff+0.3)*((W/h)+0.264)) / ((eps_reff-0.258)*((W/h)+0.8)); % extended incremental length
    
    L = (c / (2 * freq * sqrt(eps_reff))) - 2*delta_L; % length of patch

    L_eff = L + 2*delta_L; % effective antenna length

    lambda = c / freq; % wavelength

    if (W < lambda) 
        G_1 = (1/90) * (W / lambda)^2; % conductance
    else
        G_1 = (1/120) * (W / lambda)^2;
    end;

    k_0 = (2*pi) / lambda; % free space wavenumber

    integrand = @(th) ...
        (k_0*W/2)^2 .* ...
        sinc((k_0*W/2).*cos(th)/pi).^2 .* ...
        besselj(0, k_0*L.*sin(th)) .* sin(th).^3;

    G_12 = (1/(120*pi^2)) * ...
        integral(integrand, 0, pi, ...
        'RelTol', 1e-8, 'AbsTol', 1e-12); % mutual conductance

    R_in_y = 1 / (2*(G_1 + G_12)); % resonant input resistance

    y_0 = (L/pi) * acos(sqrt(R_in_y0 / R_in_y)); % inset feed point distance

    H_prime = ((Z_0 * sqrt(2 * (eps_r+1)))/119.9) + 0.5 * ((eps_r-1)/(eps_r+1)) * (log(pi/2) + (1/eps_r)*log(4/pi));
    w_h_ratio = ((exp(H_prime) / 8) - (1 / (4*exp(H_prime))))^-1; % width to height ratio of feed line

    W_0 = w_h_ratio * h; % width of the feed line

    SL = L + 6*h; % substrate/ground length and width 
    SW = W + 6*h;

    % Display results in mm
    fprintf('Substrate height: %.1f mm\n', h*1e3);
    fprintf('Substrate Length: %.3f mm\n', SL*1e3);
    fprintf('Substrate Width: %.3f mm\n', SW*1e3);
    fprintf('Effective permittivity: %.4f\n', eps_reff);
    fprintf('Patch length L: %.3f mm\n', L*1e3);
    fprintf('Effective length: %.3f mm\n', L_eff*1e3);
    fprintf('Length extension: %.3f mm\n', delta_L*1e3);
    fprintf('G1: %.6e S\n', G_1);
    fprintf('G12: %.6e S\n', G_12);
    fprintf('Edge input resistance: %.2f ohms\n', R_in_y);

    if isnan(y_0)
        fprintf('Inset feed position: no real solution\n');
    else
        fprintf('Inset feed position y0: %.3f mm\n', y_0*1e3);
    end

    fprintf('Feed-line width: %.3f mm\n\n', W_0*1e3);
end;