freq = 4.9e9;
c = 3e8;
eps_r = 4.3;

W = (c / (2*freq)) * sqrt(2 / (eps_r+1))

h_values = [0.8, 1.6, 2.4] * 1e-3;

for i = 1:length(h_values)
    h = h_values(i);

    eps_reff = ((eps_r+1) / 2) + ((eps_r-1)/2)*(1 + 12*(h/W))^(-1/2)

    delta_L = 1e3 * 0.412 * h * ((eps_reff+0.3)*((W/h)+0.264)) / ((eps_reff-0.258)*((W/h)+0.8))
    
    L = 1e3 * (c / (2 * freq * sqrt(eps_reff))) - 2*delta_L

    L_eff = L + 2*delta_L
end;