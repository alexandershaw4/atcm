function pong_tcm_active_inference_demo
% PONG_TCM_ACTIVE_INFERENCE_DEMO
%
% Minimal in-silico Pong agent whose behaviour is modulated by a
% thalamo–cortical neural mass model.
%
% Now in *classic Pong* geometry:
%   - Ball moves in 2D.
%   - Paddle is a vertical bar on the LEFT wall, moving up/down.
%
% State:
%   x_env = [ball_x; ball_y; ball_vx; ball_vy; paddle_y];

% --- Simulation parameters ------------------------------------------------
T_steps        = 500;     % number of discrete Pong time steps
dt_env         = 0.02;    % environment time step (seconds)
dt_tcm         = 1/600;    % TCM integration step (seconds)
n_tcm_substeps = 10;      % TCM micro-steps per environment step

% Pong world parameters
env.width        = 1.0;    % horizontal extent (0..1)
env.height       = 1.0;    % vertical extent (0..1)

env.paddle_x     = 0.1;    % FIXED x-position of paddle (left wall)
env.paddle_h     = 0.15;   % paddle half-height (vertical size)
env.paddle_speed = 0.6;    % paddle speed in y (units / second)

env.ball_speed   = 0.6;    % base ball speed

% Actions (discrete)
A_UP   = 1;   % move paddle UP (increase y)
A_STAY = 2;   % no move
A_DOWN = 3;   % move paddle DOWN (decrease y)
actions = [A_UP, A_STAY, A_DOWN];

% --- Initialise environment state ----------------------------------------
% State vector for environment:
% x_env = [ball_x; ball_y; ball_vx; ball_vy; paddle_y];
x_env = zeros(5,1);

x_env(1) = 0.8;                        % ball_x (start near right)
x_env(2) = 0.5;                        % ball_y (middle)

theta    = pi + (pi/3)*(rand - 0.5);   % broadly leftwards
x_env(3) = env.ball_speed * cos(theta);
x_env(4) = env.ball_speed * sin(theta);

x_env(5) = 0.5;                        % paddle_y (middle)

% --- Initialise TCM state -------------------------------------------------
[P_tcm, M_tcm, x_tcm0] = init_tcm_demo();
x_tcm = x_tcm0(:);
x_tcm = x_tcm*0;
x_tcm(1:8)=-50;

% Motor read-out from TCM: linear projection (e.g. L2/3 + L5)
% W_motor = zeros(1, numel(x_tcm));
% W_motor(2) = 0.2;
% W_motor(4) = 1.0;

% Motor read-out from TCM: start small, will be learned online
% W_motor = 1e-3 * randn(1, numel(x_tcm));
% eta_W   = 1e-3;   % learning rate for motor weights

Nx    = numel(x_tcm);

% Indices you want to use as motor readout (e.g. L2/3 + L5 pyramids)
idx_motor = [2 4 10 12];      % <<< replace with your actual indices

% Initialise W_motor only on those states
W_motor = zeros(1, Nx);
W_motor(idx_motor) = 1e-3 * randn(1, numel(idx_motor));

eta_W   = 1e-3;         % learning rate

% --- parameter vector for TCM meta-learning ---
theta   = pack_theta_from_P(P_tcm);
eta_th  = 1e-4;        % learning rate for theta (small!)
eps_th  = 1e-4;        % finite-difference step size


% --- Logging -------------------------------------------------------------
traj.x_env  = zeros(5, T_steps);
traj.x_tcm1 = zeros(1, T_steps);   % motor drive for plotting
traj.action = zeros(1, T_steps);

% --- Main simulation loop -------------------------------------------------
a_t = A_STAY;   % initial action

for t = 1:T_steps

    % 1) Update true environment given current action
    x_env = step_pong_env(x_env, a_t, env, dt_env);

    % 2) Live animation
    %draw_pong_frame(x_env, env);
    %draw_pong_frame(x_env, env, v_pred, v_teacher, a_t, theta, t*dt_env);

    % 3) Build sensory summary for TCM from game state
    s_t = build_sensory_summary(x_env);   % s_t(1) = dy = ball_y - paddle_y

    % 4) Drive TCM with this sensory input (integrate a short window)
    x_tcm = step_tcm(x_tcm, s_t, P_tcm, M_tcm, dt_tcm, n_tcm_substeps);

    % 5) "Teacher" motor drive from simple dy controller
    dy        = s_t(1);
    v_teacher = tanh(3 * dy);        % the thing you know works

    % 6) TCM-based motor drive and online learning of W_motor
    % v_pred = tanh(W_motor * x_tcm);  % current TCM->motor mapping
    % err    = v_teacher - v_pred;     % scalar error
    % 
    % % Delta rule update
    % W_motor = W_motor + eta_W * err * x_tcm';  % 1×N + (1×1 * 1×N)

    v_pred = tanh(W_motor * x_tcm);   % 1×N * N×1 -> scalar
    err    = v_teacher - v_pred;

     % --- Learn W_motor on chosen states (your sparse update) -----------
    x_m      = x_tcm(idx_motor);
    W_update = eta_W * err * x_m';           % 1×k
    W_motor(idx_motor) = W_motor(idx_motor) + W_update;

    % --- Meta-learning: update selected TCM params theta ---------------
    % Only do this every few steps if you like, to save compute:
    if mod(t, 5) == 0          % e.g. update every 5 env steps
        % Pack current theta from P_tcm (in case we've changed P elsewhere)
        theta = pack_theta_from_P(P_tcm);

        n_th   = numel(theta);
        d_v_dtheta = zeros(n_th,1);

        % compute gradient of v_pred wrt each theta_k
        for k = 1:n_th

            theta_plus      = theta;
            theta_plus(k)   = theta_plus(k) + exp(eps_th);

            % build a perturbed P
            P_plus = unpack_theta_into_P(P_tcm, theta_plus);

            % simulate TCM forward *one* micro-step with P_plus
            x_tmp = x_tcm;    % copy current state
            u     = zeros(M_tcm.ninputs,1);
            kk    = min(numel(s_t), M_tcm.ninputs);
            u(1:kk) = s_t(1:kk);

            %for it = 1:100
             
                for i = 1:n_tcm_substeps
                [dx_plus,~,~] = atcm.tc_hilge2(x_tmp, u(1), P_plus, M_tcm);
                dx_plus = real(dx_plus);
                x_tmp        = x_tmp + dt_tcm * dx_plus;
                end
            %end

            % motor drive under perturbed param
            v_plus       = tanh(W_motor * x_tmp);

            d_v_dtheta(k) = (v_plus - v_pred) / eps_th;
        end

        % Gradient descent step to reduce motor error
        % theta_new = theta + eta_th * err * d_v_dtheta
        %theta = theta + eta_th * err * real(d_v_dtheta);
        theta = theta + eta_th * -d_v_dtheta;

        theta = real(theta);

        % Write back into P_tcm
        P_tcm = unpack_theta_into_P(P_tcm, theta);
    end

    % % Only update the chosen motor states
    % x_m        = x_tcm(idx_motor);          % column
    % W_update   = eta_W * err * x_m';        % 1×k
    % W_motor(idx_motor) = W_motor(idx_motor) + W_update;

    v_des = v_pred;                  % use TCM-based drive for control

    % 7) Select next action using active-inference-flavoured scheme
    a_t = select_action_pong_tcm(x_env, v_des, env, dt_env, actions);

    draw_pong_frame(x_env, env, v_pred, v_teacher, a_t, theta, t*dt_env);

    % 8) Log
    traj.x_env(:,t)  = x_env;
    traj.x_tcm1(t)   = v_des;        % log TCM drive
    traj.action(t)   = a_t;
end

% --- Plot results --------------------------------------------------------
plot_results_pong_tcm(traj, env, dt_env);

end


function draw_pong_frame(x_env, env, v_pred, v_teacher, a_t, theta, t)
% DRAW_PONG_FRAME_EXTENDED
% Multi-panel live dashboard:
%   1. Live Pong animation
%   2. Motor drive (TCM)
%   3. Teacher vs predicted motor drive
%   4. Action trace
%   5. Parameter trajectories

persistent hFig hAx

% Extract environment info
bx = x_env(1);
by = x_env(2);
py = x_env(5);
px = env.paddle_x;

% Initialise figure on first call
if isempty(hFig) || ~isvalid(hFig)

    hFig = figure('Color','w');
    set(hFig,'Position',[200 100 900 800]);

    % Subplot layout (5 rows)
    hAx.game     = subplot(3,4,[1 2 5 6]);
    hAx.motor    = subplot(3,4,[3 4]);
    hAx.teacher  = subplot(3,4,[7 8]);
    hAx.action   = subplot(3,4,[9 10]);
    hAx.params   = subplot(3,4,[11 12]);

    % 1. GAME PANEL --------------------------------------------------------
    subplot(hAx.game); hold on;
    rectangle('Position',[0 0 env.width env.height], ...
              'EdgeColor',[0.6 0.6 0.6], 'LineWidth',1.5);

    hAx.ball = plot(bx, by, 'ko', 'MarkerFaceColor','k', 'MarkerSize',10);
    hAx.paddle = plot([px px], [py-env.paddle_h, py+env.paddle_h], ...
                      'r-', 'LineWidth',6);

    axis([0 env.width 0 env.height]); axis manual;
    set(gca,'YDir','normal'); xlabel('x'); ylabel('y');
    title('Live Pong-TCM Animation');

    % 2. MOTOR PANEL -------------------------------------------------------
    subplot(hAx.motor); hold on;
    hAx.motor_line = plot(nan, nan, 'b', 'LineWidth',1.3);
    ylabel('v_{TCM}');
    title('TCM Motor Drive');

    % 3. TEACHER vs PRED MOTOR DRIVE --------------------------------------
    subplot(hAx.teacher); hold on;
    hAx.teacher_line = plot(nan, nan, 'k', 'LineWidth',0.8);
    hAx.pred_line    = plot(nan, nan, 'r', 'LineWidth',1.2);
    ylabel('drive');
    legend({'teacher','tcm'}, 'Location','eastoutside');
    title('Teacher vs TCM Motor Drive');

    % 4. ACTION PANEL ------------------------------------------------------
    subplot(hAx.action); hold on;
    hAx.action_line = plot(nan, nan, 'b', 'LineWidth',1.2);
    yticks([1 2 3]); yticklabels({'UP','STAY','DOWN'});
    ylabel('action');
    title('Actions');

    % 5. PARAMETER PANEL ---------------------------------------------------
    subplot(hAx.params); hold on;

    % number of parameters in current theta
    n_th = numel(theta);
    hAx.n_theta = n_th;

    % colour set for different params
    cols = lines(n_th);

    % preallocate line handles
    hAx.theta_lines = gobjects(n_th,1);

    for k = 1:n_th
        hAx.theta_lines(k) = plot(nan, nan, 'LineWidth',1.2, ...
            'Color', cols(k,:));
    end

    ylabel('\theta');
    % generic legend labels: θ₁, θ₂, ...
    legStr = arrayfun(@(k) sprintf('\\theta_%d',k), 1:n_th, 'UniformOutput', false);
    legend(hAx.theta_lines, legStr, 'Location','eastoutside');
    title('TCM Parameter Learning');

    % Store handles for reuse
    guidata(hFig, hAx);
end

% Retrieve handle state
hAx = guidata(hFig);

% -------------------------------------------------------------------------
% UPDATE EACH PANEL
% -------------------------------------------------------------------------

% Update Pong animation
set(hAx.ball,   'XData', bx, ...
                'YData', by);

set(hAx.paddle, 'XData', [px px], ...
                'YData', [py-env.paddle_h, py+env.paddle_h]);


% Update motor drive trace
subplot(hAx.motor);
x = get(hAx.motor_line,'XData'); y = get(hAx.motor_line,'YData');
set(hAx.motor_line,'XData',[x t],'YData',real([y v_pred]));


% Update teacher vs TCM motor drive
subplot(hAx.teacher);
x1 = get(hAx.teacher_line,'XData'); y1 = get(hAx.teacher_line,'YData');
x2 = get(hAx.pred_line,'XData');    y2 = get(hAx.pred_line,'YData');
set(hAx.teacher_line,'XData',[x1 t],'YData',real([y1 v_teacher]));
set(hAx.pred_line,   'XData',[x2 t],'YData',real([y2 v_pred]));


% Update action trace
subplot(hAx.action);
xa = get(hAx.action_line,'XData'); ya = get(hAx.action_line,'YData');
set(hAx.action_line,'XData',[xa t],'YData',real([ya a_t]));


% Update learned parameters
subplot(hAx.params);
n_th = hAx.n_theta;

for k = 1:n_th
    hL = hAx.theta_lines(k);
    xk = get(hL,'XData');
    yk = get(hL,'YData');
    set(hL,'XData',[xk t], 'YData',real([yk theta(k)]));
end
% subplot(hAx.params);
% xp1 = get(hAx.th1,'XData'); yp1 = get(hAx.th1,'YData');
% xp2 = get(hAx.th2,'XData'); yp2 = get(hAx.th2,'YData');
% 
% set(hAx.th1,'XData',[xp1 t],'YData',[yp1 theta(1)]);
% set(hAx.th2,'XData',[xp2 t],'YData',[yp2 theta(2)]);

drawnow limitrate;
end



% =====================================================================
%              ENVIRONMENT + DYNAMICS FUNCTIONS
% =====================================================================
function theta = pack_theta_from_P(P)
% PACK_THETA_FROM_P  Pick out the TCM parameters we want to learn.
%
% Example: a scalar input gain and a single connectivity weight.
% Adjust this to use the fields you really care about.


theta(1) = exp(P.H(2,2));        % <- replace with your actual field
theta(2) = exp(P.H(4,4));      % <- and this

theta(3) = exp(P.Hn(2,2));        % <- replace with your actual field
theta(4) = exp(P.Hn(4,4));      % <- and this

theta(5) = exp(P.H(1,8));        % <- replace with your actual field


end

function P = unpack_theta_into_P(P, theta)
% UNPACK_THETA_INTO_P  Write learned params back into P.

P.H(2,2) = log(theta(1));      % same mapping as above
P.H(4,4) = log(theta(2));

P.Hn(2,2) = log(theta(3));      % same mapping as above
P.Hn(4,4) = log(theta(4));

P.H(1,8) = log(theta(5));      % <- replace with your actual field


end



function x_env_next = step_pong_env(x_env, a_t, env, dt)
% STEP_PONG_ENV  One-step update of classic Pong world.
%
% x_env = [ball_x; ball_y; ball_vx; ball_vy; paddle_y];

bx  = x_env(1);
by  = x_env(2);
vx  = x_env(3);
vy  = x_env(4);
py  = x_env(5);

% --- Update paddle position (vertical movement) ---
switch a_t
    case 1 % A_UP -> move paddle up (increase y)
        py = py + env.paddle_speed * dt;
    case 2 % A_STAY
        % no change
    case 3 % A_DOWN -> move paddle down (decrease y)
        py = py - env.paddle_speed * dt;
end

% enforce vertical boundaries (0..height)
py = max(0, min(env.height, py));

% --- Update ball position ---
bx = bx + vx * dt;
by = by + vy * dt;

% Reflect ball on top/bottom walls
if by <= 0
    by = -by;
    vy = -vy;
elseif by >= env.height
    by = 2*env.height - by;
    vy = -vy;
end

% Reflect ball on right wall
if bx >= env.width
    bx = 2*env.width - bx;
    vx = -vx;
end

% --- Check collision with paddle (vertical bar at x = paddle_x) ---
if bx <= env.paddle_x
    % Paddle spans [py - paddle_h, py + paddle_h] in y
    if by >= py - env.paddle_h && by <= py + env.paddle_h
        % Bounce: mirror around paddle_x and flip vx
        bx = 2*env.paddle_x - bx;
        vx = -vx;
    else
        % Miss: reset ball from right with random direction
        bx = 0.9 * env.width;
        by = rand * env.height;
        theta = pi + (pi/2)*(rand - 0.5);  % broadly leftwards
        vx = env.ball_speed * cos(theta);
        vy = env.ball_speed * sin(theta);
    end
end

x_env_next = [bx; by; vx; vy; py];

end


function s_t = build_sensory_summary(x_env)
% BUILD_SENSORY_SUMMARY  Low-dimensional "sensory" projection for TCM drive.
%
% Here we use:
%   dy   = ball_y - paddle_y  (relevant error)
%   xb   = ball_x
%   vby  = ball_vy

bx = x_env(1);
by = x_env(2);
vy = x_env(4);
py = x_env(5);

dy = by - py;

s_t = [dy; bx; vy];

end


% =====================================================================
%                       TCM INTEGRATION
% =====================================================================
function x_tcm = step_tcm(x_tcm, s_t, P, M, dt, n_sub)
% STEP_TCM  Integrate thalamo–cortical model for a few micro-steps.
%
% Uses:
%   [dx,J,D] = atcm.tc_hilge2(x,u,P,M)

n_inputs = M.ninputs;
u        = zeros(n_inputs,1);
k        = min(numel(s_t), n_inputs);
u(1:k)   = s_t(1:k);

x_tcm = x_tcm(:);

for i = 1:n_sub
    [dx,~,~] = atcm.tc_hilge2(x_tcm, u, P, M);
    x_tcm    = x_tcm + dt * dx;
end

end





function [P, M, x] = init_tcm_demo()
% INIT_TCM_DEMO  Load TCM parameters from tcm_init.mat
p = mfilename('fullpath');
[p,~,~] = fileparts(p);

load([p '/tcm_init.mat']);      % must contain tcm_init.pE and tcm_init.M

P = tcm_init.pE;
M = tcm_init.M;
x = M.x;

% Make sure ninputs is set appropriately
if ~isfield(M,'ninputs') || isempty(M.ninputs)
    M.ninputs = 3;   % to match s_t = [dy; bx; vy]
end

end


% =====================================================================
%                       ACTION SELECTION
% =====================================================================

% function draw_pong_frame(x_env, env)
% % DRAW_PONG_FRAME  Live 2D animation of Pong world (classic Pong).
% %
% % x_env = [ball_x; ball_y; ball_vx; ball_vy; paddle_y];
% 
% persistent hFig hBall hPaddle isInit
% 
% bx = x_env(1);
% by = x_env(2);
% py = x_env(5);
% px = env.paddle_x;
% 
% % --- INITIALISE FIGURE ON FIRST CALL ---
% if isempty(isInit) || ~isvalid(hFig)
%     hFig = figure('Color','w');
%     set(hFig,'Position',[100 100 400 400]);
%     hold on;
% 
%     % Court boundary
%     rectangle('Position',[0 0 env.width env.height], ...
%               'EdgeColor',[0.6 0.6 0.6], 'LineWidth',1.5);
% 
%     % Ball marker
%     hBall = plot(bx, by, 'ko', 'MarkerFaceColor','k', ...
%                      'MarkerSize',10);
% 
%     % Paddle as vertical bar
%     hPaddle = plot([px px], [py-env.paddle_h, py+env.paddle_h], ...
%                    'r-', 'LineWidth',6);
% 
%     axis([0 env.width 0 env.height]);
%     axis manual;
%     set(gca,'YDir','normal');
%     xlabel('x'); ylabel('y');
%     title('Live Pong-TCM Animation');
% 
%     isInit = true;
% else
%     % --- UPDATE MARKERS ---
%     set(hBall,   'XData', bx, 'YData', by);
%     set(hPaddle, 'XData', [px px], ...
%                  'YData', [py-env.paddle_h, py+env.paddle_h]);
% end
% 
% drawnow limitrate;   % keeps animation responsive without lagging CPU
% end
% 
function a_opt = select_action_pong_tcm(x_env, v_des, env, dt_env, actions)

lambda_motor = 0.5;

costs = zeros(numel(actions),1);

for i = 1:numel(actions)
    a = actions(i);

    x_next = step_pong_env(x_env, a, env, dt_env);

    by_next = x_next(2);
    py_next = x_next(5);

    align_cost = (by_next - py_next)^2;

    switch a
        case 1, u_dir = +1;  % UP
        case 2, u_dir = 0;   % STAY
        case 3, u_dir = -1;  % DOWN
    end

    desired_dir = sign(v_des);
    motor_cost  = (u_dir - desired_dir).^2;

    costs(i) = align_cost + lambda_motor * motor_cost;
end

[~, idx] = min(costs);
a_opt    = actions(idx);
end




% =====================================================================
%                           PLOTTING
% =====================================================================

function plot_results_pong_tcm(traj, env, dt_env)

t = (0:size(traj.x_env,2)-1) * dt_env;

bx = traj.x_env(1,:);
by = traj.x_env(2,:);
py = traj.x_env(5,:);
v  = traj.x_tcm1;

figure;
subplot(3,1,1);
plot(t, by, 'LineWidth',1.5); hold on;
plot(t, py, 'LineWidth',1.5);
xlabel('t');
ylabel('y');
legend({'ball y','paddle y'});
title('Vertical positions');

subplot(3,1,2);
plot(t, v, 'LineWidth',1.5);
xlabel('t');
ylabel('motor drive');
title('TCM-derived motor signal');

subplot(3,1,3);
plot(t, traj.action, 'LineWidth',1.5);
xlabel('t');
ylabel('action');
title('Actions over time');

end
