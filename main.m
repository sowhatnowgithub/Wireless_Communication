%% ======================================================================
% Wireless Communication Project.
%
% Title : UAV-Enabled Secure ISAC 
% Name : J L S Pavan Kumar
% Roll No : 2302cm05
%% ======================================================================
clear; clc; close all;

%% PARAMETERS
M = 4;
D = 200;
T = 12;
N = 24;
ts = T / N;
fc = 6e9;
lambda_c = 3e8 / fc;
d_ant = lambda_c / 2;

beta0 = 10^(-30/10)*1e-3;
sigma2 = 10^(-90/10)*1e-3;
Pmax = 10^(30/10)*1e-3;

vmax = 25;
Vmax = vmax * ts;

Gamma_t = 1e-6;
Gamma_e = 1e-6;

%% NODE LOCATIONS
s_u = [250; 520];
s_t = [250; 480];
s_e = [350; 500];
rho_I = [300; 400];
rho_F = [300; 600];

%% ITERATION SETTINGS
max_AO_iter = 3;
max_SCA1_iter = 8;
step_size_traj = 5;   
max_traj_iter = 30;


%% HELPER FUNCTIONS 
steer_fn = @(M_in, rho_uav, s_gnd) ...
    exp(1j*2*pi*d_ant/lambda_c*(0:M_in-1).' * ...
    (norm(rho_uav-s_gnd)/sqrt(norm(rho_uav-s_gnd)^2+D^2)));

chan_fn = @(M_in, rho_uav, s_gnd) ...
    sqrt(beta0/(norm(rho_uav-s_gnd)^2+D^2)) * steer_fn(M_in, rho_uav, s_gnd);

% Curved initial trajectory: arc around the eavesdropper
rho = zeros(2,N);
for n = 1:N
    t = (n-1)/(N-1);
    rho(1,n) = rho_I(1) + t*(rho_F(1)-rho_I(1))  - 80*sin(pi*t);
    rho(2,n) = rho_I(2) + t*(rho_F(2)-rho_I(2)) + 30*sin(pi*t);
end


%% INITIALIZE BEAMFORMING
B_sol = cell(1,N); As_sol = cell(1,N);
for n = 1:N
    B_sol{n}  = (Pmax/(2*M))*eye(M);
    As_sol{n} = (Pmax/(2*M))*eye(M);
end

objective_history = zeros(max_AO_iter, 1);
fprintf('===== Starting AO =====\n');

for iter_AO = 1:max_AO_iter
    fprintf('\n--- AO Iteration %d/%d ---\n', iter_AO, max_AO_iter);
    
    %% BEAMFORMING (SDR + SCA) 
    B_local = B_sol; As_local = As_sol;
    
    for iter_s1 = 1:max_SCA1_iter
        for n = 1:N
            h_u = chan_fn(M, rho(:,n), s_u);
            h_e = chan_fn(M, rho(:,n), s_e);
            a_t = steer_fn(M, rho(:,n), s_t);
            a_e = steer_fn(M, rho(:,n), s_e);
            dt2 = norm(rho(:,n)-s_t)^2 + D^2;
            de2 = norm(rho(:,n)-s_e)^2 + D^2;
            Hu = h_u*h_u'; He = h_e*h_e';
           
            nu_l = real(trace(Hu*B_local{n}));
            du_l = real(trace(Hu*As_local{n})) + sigma2;
            ne_l = real(trace(He*B_local{n}));
            de_l = real(trace(He*As_local{n})) + sigma2;
            
            cu1 = 1/(log(2)*(nu_l+du_l));
            cu2 = nu_l/(log(2)*du_l*(nu_l+du_l));
            ce1 = 1/(log(2)*(ne_l+de_l));
            ce2 = ne_l/(log(2)*de_l*(ne_l+de_l));
            
            Rc_u = log2(1+nu_l/du_l) - cu1*nu_l + cu2*real(trace(Hu*As_local{n}));
            Rc_e = log2(1+ne_l/de_l) - ce1*ne_l + ce2*real(trace(He*As_local{n}));
            
            cvx_begin sdp quiet
                cvx_solver mosek
                variable Bv(M,M) hermitian semidefinite
                variable Asv(M,M) hermitian semidefinite
                Ru_s = Rc_u + cu1*real(trace(Hu*Bv)) - cu2*real(trace(Hu*Asv));
                Re_s = Rc_e + ce1*real(trace(He*Bv)) - ce2*real(trace(He*Asv));
                
                maximize(Ru_s - Re_s  )

                subject to
                    real(trace(Bv))+real(trace(Asv)) <= Pmax;
                    real(a_t'*(Bv+Asv)*a_t) >= Gamma_t*dt2;
                    real(a_e'*(Bv+Asv)*a_e) <= Gamma_e*de2;
            cvx_end
            
            if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
                B_sol{n} = Bv; As_sol{n} = Asv;
            end
        end
        B_local = B_sol; As_local = As_sol;
    end
    
    % Rank-1 recovery
    for n = 1:N
        [V,Dg] = eig(B_sol{n}); dv = real(diag(Dg));
        [mv,idx] = max(dv);
        if mv > 0; b_n = sqrt(mv)*V(:,idx); else; b_n = zeros(M,1); end
        B_sol{n} = b_n*b_n';
    end
    
    fprintf('  Beamforming done.\n');
    
    %%  TRAJECTORY (GRADIENT PROJECTION) 
    step_sz = step_size_traj;
    
    for iter_t = 1:max_traj_iter
        grad_rho = zeros(2, N);
        R_current = compute_secrecy(rho, B_sol, As_sol, N, M, ...
            s_u, s_e, sigma2, beta0, d_ant, lambda_c, D);
        
        for n = 2:N-1
            for dim = 1:2
                rho_p = rho; rho_m = rho;
                rho_p(dim,n) = rho(dim,n) + 0.3;
                rho_m(dim,n) = rho(dim,n) - 0.3;
                
                Rp = compute_secrecy(rho_p, B_sol, As_sol, N, M, ...
                    s_u, s_e, sigma2, beta0, d_ant, lambda_c, D);
                Rm = compute_secrecy(rho_m, B_sol, As_sol, N, M, ...
                    s_u, s_e, sigma2, beta0, d_ant, lambda_c, D);
                
                grad_rho(dim,n) =(Rp - Rm) / 1.0;
            end
        end
        
        % Gradient step
        rho_new = rho + step_sz * grad_rho;
        rho_new(:,1) = rho_I;
        rho_new(:,N) = rho_F;
        
        % Project: enforce speed constraints (forward pass)
        for n = 2:N
            disp_vec = rho_new(:,n) - rho_new(:,n-1);
            if norm(disp_vec) > Vmax
                rho_new(:,n) = rho_new(:,n-1) + Vmax*disp_vec/norm(disp_vec);
            end
        end
        % Backward pass
        for n = N-1:-1:1
            disp_vec = rho_new(:,n) - rho_new(:,n+1);
            if norm(disp_vec) > Vmax
                rho_new(:,n) = rho_new(:,n+1) + Vmax*disp_vec/norm(disp_vec);
            end
        end
        rho_new(:,1) = rho_I;
        rho_new(:,N) = rho_F;
        
        % Check improvement
        R_new = compute_secrecy(rho_new, B_sol, As_sol, N, M, ...
            s_u, s_e, sigma2, beta0, d_ant, lambda_c, D);
        
        if R_new >= R_current
            rho = rho_new;
            fprintf('  Traj iter %d: %.4f -> %.4f\n', iter_t, R_current, R_new);
        else
            step_sz = step_sz * 0.5;
            if step_sz < 0.01
                fprintf('  Trajectory converged at iter %d\n', iter_t);
                break;
            end
        end
    end
    
    %% COMPUTE OBJECTIVE
    avg_secrecy = compute_secrecy(rho, B_sol, As_sol, N, M, ...
        s_u, s_e, sigma2, beta0, d_ant, lambda_c, D);
    objective_history(iter_AO) = avg_secrecy;
    fprintf('  => Avg Secrecy Rate = %.4f bps/Hz\n', avg_secrecy);
    
    if iter_AO > 1 && abs(objective_history(iter_AO)-objective_history(iter_AO-1)) < 1e-3
        fprintf('\nAO Converged at iteration %d\n', iter_AO);
        break;
    end
end

%%  PLOTS 
rho_proposed = rho;
B_proposed = B_sol;
As_proposed = As_sol;
%% Fig 0: Convergence
figure;
valid = objective_history > 0;
plot(find(valid), objective_history(valid), 'b-o', 'LineWidth',1.5);
xlabel('AO Iteration'); ylabel('Avg Secrecy Rate (bps/Hz)');
title('Convergence'); grid on;
%% Fig 1: Trajectory
figure;
plot(rho_proposed(1,:), rho_proposed(2,:), 'b-o', 'LineWidth',1.5,'MarkerSize',6); hold on;
plot(s_u(1),s_u(2),'gs','MarkerSize',14,'MarkerFaceColor','g');
plot(s_t(1),s_t(2),'b^','MarkerSize',14,'MarkerFaceColor','b');
plot(s_e(1),s_e(2),'rx','MarkerSize',14,'LineWidth',3);
plot(rho_I(1),rho_I(2),'kp','MarkerSize',18,'MarkerFaceColor','k');
plot(rho_F(1),rho_F(2),'kd','MarkerSize',18,'MarkerFaceColor','k');
legend('Trajectory','User','Target','Eavesdropper','Start','End');
xlabel('x (m)'); ylabel('y (m)');
title('Fig 1: UAV Trajectory'); grid on; axis equal; hold off;

%% Fig 2: Power Allocation 

power_c = zeros(1,N); 
power_s = zeros(1,N);
power_total = zeros(1,N);

for n = 1:N
    power_c(n) = real(trace(B_proposed{n}));
    power_s(n) = real(trace(As_proposed{n}));
    power_total(n) = power_c(n) + power_s(n);
end

figure;

plot(1:N, power_c, 'b-o', 'LineWidth',2, 'MarkerSize',6); hold on;
plot(1:N, power_s, 'r-s', 'LineWidth',2, 'MarkerSize',6);
plot(1:N, power_total, 'k--', 'LineWidth',2);

xlabel('Time Slot');
ylabel('Power (W)');
title('Fig 2: Power Allocation (Communication vs Sensing)');
legend('Communication Power','Sensing (AN) Power','Total Power','Location','best');

grid on;
ylim([0, Pmax*1.1]);   % ensures starting from 0 clearly

%%  Fig 3: Beampattern

n_plot = min(5, N);

% GRID
x_g = linspace(0,600,80); 
y_g = linspace(0,600,80);
[Xg, Yg] = meshgrid(x_g, y_g);
gain = zeros(size(Xg));

% BEAMPATTERN 
for ii = 1:length(y_g)
    for jj = 1:length(x_g)
        a_test = steer_fn(M, rho_proposed(:,n_plot), [Xg(ii,jj);Yg(ii,jj)]);
        gain(ii,jj) = real(a_test'*(B_proposed{n_plot}+As_proposed{n_plot})*a_test);
    end
end

% METRICS 

% Channels (for SINR)
h_u = chan_fn(M, rho_proposed(:,n_plot), s_u);
h_e = chan_fn(M, rho_proposed(:,n_plot), s_e);

% Steering (for gains)
a_u = steer_fn(M, rho_proposed(:,n_plot), s_u);
a_t = steer_fn(M, rho_proposed(:,n_plot), s_t);

Bv = B_proposed{n_plot};
Asv = As_proposed{n_plot};

% --- USER ---
gain_u = real(a_u' * Bv * a_u);   % FIXED (array gain)
sinr_u = real(h_u' * Bv * h_u) / (real(h_u' * Asv * h_u) + sigma2);

% --- EAVESDROPPER ---
sinr_e = real(h_e' * Bv * h_e) / (real(h_e' * Asv * h_e) + sigma2);

% --- TARGET ---
gain_t = real(a_t' * (Bv + Asv) * a_t);
%  PLOT 
figure;

% Beampattern
surf(Xg, Yg, 10*log10(gain + 1e-20), 'EdgeColor', 'none');
view(2);
colorbar;
colormap('jet');
hold on;

% Plot nodes
plot3(s_u(1), s_u(2), 50, 'gs', 'MarkerSize',14, 'MarkerFaceColor','g');
plot3(s_t(1), s_t(2), 50, 'b^', 'MarkerSize',14, 'MarkerFaceColor','b');
plot3(s_e(1), s_e(2), 50, 'rx', 'MarkerSize',14, 'LineWidth',3);

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Fig 4: Beampattern (Slot %d)', n_plot));

% TEXT ANNOTATION 

text(s_u(1)+10, s_u(2)+10, 60, ...
    sprintf('User\nGain = %.3f\nSINR = %.3f', gain_u, sinr_u), ...
    'Color','b','FontSize',10,'FontWeight','bold');
% EAVESDROPPER
text(s_e(1)+10, s_e(2)+10, 60, ...
    sprintf('Eve\nSINR = %.3f', sinr_e), ...
    'Color','b','FontSize',10,'FontWeight','bold');

% TARGET
text(s_t(1)+10, s_t(2)+10, 60, ...
    sprintf('Target\nGain = %.3f', gain_t), ...
    'Color','b','FontSize',10,'FontWeight','bold');

hold off;

%% Fig 4: Secrecy vs M

fprintf('\n===== Monte Carlo for Fig 3 =====\n');

M_vals = [4, 6, 8, 10, 12, 14, 16, 18];    
num_MC = 5;                

avg_sec_proposed = zeros(length(M_vals),1);
avg_sec_straight = zeros(length(M_vals),1);

for mi = 1:length(M_vals)

    Mm = M_vals(mi);
    sec_vals_prop = zeros(num_MC,1);
    sec_vals_str  = zeros(num_MC,1);

    fprintf('M = %d: ', Mm);

    for mc = 1:num_MC

        % Random positions (same region)
        su_r = [200+200*rand; 400+200*rand];
        st_r = [200+200*rand; 400+200*rand];
        se_r = [200+200*rand; 400+200*rand];

        % Proposed (optimized trajectory)
        sec_vals_prop(mc) = run_AO_scheme(Mm,D,N,Pmax,...
            sigma2,beta0,d_ant,lambda_c,Vmax,Gamma_t,Gamma_e,...
            su_r,st_r,se_r,rho_I,rho_F,'proposed');

        % Straight trajectory baseline
        sec_vals_str(mc) = run_AO_scheme(Mm,D,N,Pmax,...
            sigma2,beta0,d_ant,lambda_c,Vmax,Gamma_t,Gamma_e,...
            su_r,st_r,se_r,rho_I,rho_F,'straight');

        fprintf('.');
    end

    avg_sec_proposed(mi) = mean(sec_vals_prop);
    avg_sec_straight(mi) = mean(sec_vals_str);

    fprintf(' done\n');
end

%% Plot 
figure;
plot(M_vals, avg_sec_proposed, 'b-o', 'LineWidth',2,'MarkerSize',8); hold on;
plot(M_vals, avg_sec_straight, 'r-^', 'LineWidth',2,'MarkerSize',8);

xlabel('Number of Antennas (M)');
ylabel('Average Secrecy Rate (bps/Hz)');
title('Fig 3: Average Secrecy Rate vs Number of Antennas');
legend('Proposed (Optimized Trajectory)','Straight Trajectory','Location','northwest');
grid on;
