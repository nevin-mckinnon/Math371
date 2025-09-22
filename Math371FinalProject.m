% Define the constants and initial conditions
alpha = 0.1; % radiobiological parameter (example value), Consider reducing slightly if needed
beta = 0.05; % the parameter for the linear-quadratic model (example value), Increase beta to emphasize the quadratic component of radiation effect

lambda_h = 0.05; % Initial growth rate for healthy cells, consider it lower than lambda_0
gamma_h = 0.005; % Reduced effectiveness of radiation on healthy cells

VH_initial = 500;
VT = 1;
K = 500;

theta = 0.015; % vascular growth retardation factor (example value)
lambda_0 = 1; % initial tumor growth rate (example value)
eta_cl = 0.15; % cell clearance rate (example value)
D = 5; % dose of radiation (example value)
tR = 47; % time when radiation is given (start of the treatment)
trad = 100; % duration of radiation effect
t_end = trad + trad + tR;
t_star = 0.1; %effect of radiation on tumour growth
comp_T_on_H = 0.8; %comp of tumour on healthy
comp_H_on_T = 3; %comp of healthy on tumour


%define the time spans and initial conditions for each phase
tspan_before = [0, tR];
% Adjusted initial conditions before radiation to include VND
initial_conditions_before = [VT; 0; lambda_0; VH_initial]; % Include initial volume of healthy cells

disp(['VH before: ', num2str(VH_initial)])


% Solve ODEs before radiation

% ODE solver for the phase before radiation
[t_before, y_before] = ode23(@(t, y) odes_before_radiation_comp(t, y, theta, lambda_0, lambda_h, K, comp_T_on_H, comp_H_on_T), tspan_before, initial_conditions_before);

% Adjusted to include VH

% Initial conditions set prior to radiation treatment
initial_conditions_during = [y_before(end, 1);  % VT from the end of 'Before Radiation'
                             0;                  % Assuming VND resets or remains unchanged initially
                             y_before(end, 3);   % lambda from the end of 'Before Radiation'
                             y_before(end, 4)];  % VH from the end of 'Before Radiation'

disp(['VH during: ', num2str(y_before(end, 4))])



% Execute the Optimization
[optimal_params, optimal_value] = runOptimizationWithDynamicInitialConditions(initial_conditions_during, tR, t_end, t_star, theta, lambda_0, eta_cl, gamma_h, lambda_h, K, comp_T_on_H, comp_H_on_T);

% Extract and Display Optimized Parameters
D = optimal_params(1);
trad = optimal_params(2);

disp(['D: ', num2str(D),', trad: ', num2str(trad)]);

% Assuming initial conditions are set appropriately before this step
% and tR, t_end, theta, lambda_0, eta_cl, t_star, gamma_h, lambda_h, K, comp_T_on_H, and comp_H_on_T are defined elsewhere in your code


% Set the simulation time span for the during radiation phase based on optimized trad
t_span_during = [tR, tR + trad];

% Solve ODEs during radiation with updated initial conditions and OPTIMIZED parameters
[t_during, y_during] = ode23(@(t, y) odes_during_radiation_comp(t, y, D, alpha, beta, theta, lambda_0, eta_cl, t_star, trad, tR, t_end, gamma_h, lambda_h, K , comp_T_on_H, comp_H_on_T), t_span_during, initial_conditions_during);

% Define the time span for the "after radiation" phase
disp(['VH after : ', num2str(y_during(end, 4))]);

tspan_after = [tR + trad, tR + trad + trad];

%initlve ODEs after radiation with updated initial conditions including VH the results

% Extracting the state at the end of the 'during radiation' phase
VT_end_during = y_during(end, 1); % VT from the end of 'During Radiation'
VND_end_during = y_during(end, 2); % VND from the end of 'During Radiation'
lambda_end_during = y_during(end, 3); % lambda from the end of 'During Radiation'
VH_end_during = y_during(end, 4); % VH from the end of 'During Radiation'

% Setting up initial conditions for the 'after radiation' phase
initial_conditions_after = [VT_end_during;  % Carry forward VT
                            VND_end_during; % Carry forward VND (assuming it's relevant post-radiation)
                            lambda_end_during; % Carry forward lambda, though it might be adjusted if you have a dynamic lambda model
                            VH_end_during]; % Carry forward VH

disp(['VH at the start of after radiation phase: ', num2str(VH_end_during)]);


% ODE solver for the phase after radiation
[t_after, y_after] = ode23(@(t, y) odes_after_radiation_comp(t, y, theta, lambda_0, eta_cl, lambda_h,K, comp_T_on_H, comp_H_on_T), tspan_after, initial_conditions_after);

t_combined = [t_before; t_during; t_after];
y_combined = [y_before; y_during; y_after];

clf; % Clear the current figure

% Stacked plot for tumor volume, non-dividing cells, and healthy cells volume
stacked_data = [y_combined(:, 1), y_combined(:, 2), y_combined(:, 4)];
h = area(t_combined, stacked_data, 'LineStyle', 'none'); % Keep handle to the area plot
xlabel('Time in Days');
ylabel('Total Volume');
title('Total Volume by Component of treated Stage 4 Cancer');
% Add vertical lines for radiation phases
hold on;
xline(tR, '--k', 'Start of Radiation', 'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'center');
xline(tR + trad, '--k', 'End of Radiation', 'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'center');
hold off;
% Update legend to only include the area plot handles
legend(h, 'Volume of tumor', 'Volume of non-dividing cells', 'Volume of healthy cells', 'Location', 'northoutside', 'Orientation', 'horizontal');

% Annotate D and trad in the top right corner of the graph
text(1, 1, sprintf('Optimal Dosage: %.1f\nLength of Dosage: %d', D, round(trad)), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Margin', 2, 'BackgroundColor', 'white');




% Function for p(D)
function p = pD(D, t_star, trad, tR, t_end, alpha, beta)
    chi = chiD(D, alpha, beta);  % Use the chi function already defined
    Tm = t_end - (tR + trad);  % Observation time after the end of radiation effect
    p = 1 - (t_star / (3 * Tm)) * chi;
end

% Function for g(D)
function g = gD(D, trad, tR, t_end, alpha, beta)
    chi = chiD(D, alpha, beta);  % Use the chi function already defined
    Tm = t_end - (tR + trad);  % Observation time after the end of radiation effect
    g = chi / (3 * Tm);
end


% Function for the linear-quadratic model chi(D)
function chi = chiD(D, alpha, beta)
    chi = alpha * D * (1 + D / (alpha / beta));
end


% ODE function before radiation, including logistic growth for healthy cells

function dydt = odes_before_radiation_comp(t, y, theta, lambda_0, lambda_h, K, comp_T_on_H, comp_H_on_T)
    VT = y(1); % Tumor volume
    VND = y(2); % Non-dividing tumor cells
    lambda = y(3); % Tumor growth rate
    VH = y(4); % Healthy cell volume

    % Modified growth rates to include competition (Lotka-Volterra model)
    dVTdt = lambda * VT * (1 - (VT + comp_T_on_H * VH) / K); % Tumor growth rate adjusted for competition
    dVHdt = lambda_h * VH * (1 - (VH + comp_H_on_T * VT) / K); % Healthy cell growth rate adjusted for competition

    % Update the differential equations to include the adjusted growth rates
    dydt = [dVTdt; % Adjusted growth of tumor cells including competition
            0; % No change in non-dividing cells before radiation
            -theta * lambda; % Adjusted growth rate due to vascular retardation
            dVHdt]; % Adjusted logistic growth of healthy cells including competition
end


% ODE function during radiation, adjusted for logistic growth of healthy cells

function dydt = odes_during_radiation_comp(t, y, D, alpha, beta, theta, lambda_0, eta_cl, t_star, trad, tR, t_end, gamma_h, lambda_h, K, comp_T_on_H, comp_H_on_T)
    VT = y(1); % Tumor volume
    VND = y(2); % Non-dividing tumor cells
    lambda = y(3); % Tumor growth rate
    VH = y(4); % Healthy cell volume

    % Radiation effects
    p = pD(D, t_star, trad, tR, t_end, alpha, beta);
    g = gD(D, trad, tR, t_end, alpha, beta);
    p_h = 1 - gamma_h * (1 - p); % Adjusted probability for healthy cells
    g_h = g*gamma_h;
    %disp(['p: ', num2str(p), ' g: ', num2str(g) ' p_h: ', num2str(p_h)]);

    % Competition model integration
    % Adjusted growth rate for tumor cells considering competition and radiation
    dVTdt = (lambda * VT * (1 - (VT + comp_T_on_H * VH) / K)) * p - g * VT;
    
    % Adjusted growth rate for healthy cells considering competition and radiation
    dVHdt = lambda_h * VH * (1 - (VH + comp_H_on_T * VT) / K) * p_h - g_h * VH;

    dydt = [dVTdt; % Adjusted growth of tumor cells including competition and radiation
            (g * VT + g_h * VH) - eta_cl * VND; % Transition to and clearance of non-dividing cells
            -theta * lambda * lambda_0; % Adjusted tumor growth rate due to vascular retardation
            dVHdt]; % Adjusted logistic growth of healthy cells including competition and radiation
end

% ODE function after radiation, adjusted for logistic growth_of healthy cells

function dydt = odes_after_radiation_comp(t, y, theta, lambda_0, eta_cl, lambda_h, K, comp_T_on_H, comp_H_on_T)
    VT = y(1); % Tumor volume
    VND = y(2); % Non-dividing tumor cells
    lambda = y(3); % Tumor growth rate
    VH = y(4); % Healthy cell volume

    % Logistic growth for healthy cells adjusted for competition
    dVHdt = lambda_h * VH * (1 - (VH + comp_H_on_T * VT) / K);
    
    % Growth rate for tumor cells adjusted for competition
    dVTdt = lambda * VT * (1 - (VT + comp_T_on_H * VH) / K);
    
    % No direct competition effect on VND, but including its dynamics for completeness
    dVNDdt = -eta_cl * VND; % Clearance of non-dividing cells

    dydt = [dVTdt; % Adjusted growth of tumor cells including competition
            dVNDdt; % Clearance of non-dividing cells
            -theta * lambda * lambda_0; % Adjusted tumor growth rate due to vascular retardation
            dVHdt]; % Adjusted logistic growth of healthy cells including competition
end

% Objective Function for Optimization


function objValue = objectiveFunction(params, initial_conditions, tR, t_end, t_star, theta, lambda_0, eta_cl, gamma_h, lambda_h, K, comp_T_on_H, comp_H_on_T)
    % Unpack parameters
    D = params(1);
    trad = params(2);

    % Constants for alpha and beta (assumed constants, not optimized)
    alpha = 0.1;
    beta = 0.05;

    % Run the simulation
    t_span = [tR, tR + trad];
    [~, y] = ode23(@(t, y) odes_during_radiation_comp(t, y, D, alpha, beta, theta, lambda_0, eta_cl, t_star, trad, tR, t_end, gamma_h, lambda_h, K, comp_T_on_H, comp_H_on_T), t_span, initial_conditions);
    
    final_VT = y(end, 1);
    initial_VT = initial_conditions(1);
    final_VH = y(end, 4);
    initial_VH = initial_conditions(4);

    % Calculate percentage change for VT and VH
    pct_change_VT = ((final_VT - initial_VT) / initial_VT) * 100;
    pct_change_VH = ((final_VH - initial_VH) / initial_VH) * 100;
    

    % Revised Objective: Maximize the decrease in VT while minimizing the decrease in VH
    % We want a high negative value for pct_change_VT and a low (or zero) negative value for pct_change_VH
    weighting_factor_VT = 0.7; % You can adjust this to make decreasing VT more significant
    weighting_factor_VH = 1; % You can adjust this to penalize decreases in VH more heavily
    
    % Adjusting the objective function to better balance the objectives

    objValue = (weighting_factor_VT * pct_change_VT) + (weighting_factor_VH * max(0, -pct_change_VH));
    % Here, we subtract the weighted pct_change_VT (to reward its decrease) and also subtract the penalty for VH decrease
    % if pct_change_VH is negative, max(0, -pct_change_VH) is positive, applying a penalty. If VH increases or remains the same, the penalty is zero.
end


% Wrapper Function for Dynamic Initial Conditions and Optimization

function [optimal_params, optimal_value] = runOptimizationWithDynamicInitialConditions(initial_conditions_during, tR, t_end, t_star, theta, lambda_0, eta_cl, gamma_h, lambda_h, K, comp_T_on_H, comp_H_on_T)
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

    objectiveFunctionWithIC = @(params) objectiveFunction(params, initial_conditions_during, tR, t_end, t_star, theta, lambda_0, eta_cl, gamma_h, lambda_h, K, comp_T_on_H, comp_H_on_T);

    initial_guess = [5, 50]; % Only D and trad
    lb = [0.5, 1]; % Lower bounds for D and trad
    ub = [100, inf]; % Upper bounds for D and trad

    [optimal_params, optimal_value] = fmincon(objectiveFunctionWithIC, initial_guess, [], [], [], [], lb, ub, [], options);
end
