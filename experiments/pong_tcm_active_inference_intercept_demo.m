function pong_tcm_active_inference_intercept_demo
% PONG_TCM_ACTIVE_INFERENCE_INTERCEPT_DEMO
%
% Drop-in replacement for the teacher-dependent Pong-TCM demo.
%
% Core change:
%   The TCM no longer learns to copy v_teacher = tanh(3*dy).
%   Instead, it learns a self-supervised intercept target:
%
%       "Where will the ball cross the paddle x-coordinate?"
%
%   The target drive is based on the difference between the paddle y-position
%   and the predicted future ball y-position at paddle contact.
%
% Architecture:
%   Pong state -> sensory summary -> TCM dynamics -> learned motor readout
%   -> short-horizon action rollout -> paddle action.
%
% Notes:
%   - TCM parameters are deliberately frozen here. This makes the TCM act as a
%     nonlinear dynamical feature generator / reservoir.
%   - Only the sparse motor readout W_motor is learned online using a
%     normalized LMS-style update.
%   - The intercept target is self-supervised from game physics, not a hand-coded
%     dy controller.
%
% Required file:
%   tcm_init.mat in the same folder, containing tcm_init.pE and tcm_init.M
%
% Required model:
%   atcm.tc_hilge2 on the MATLAB path.

% -------------------------------------------------------------------------
% Simulation parameters
% -------------------------------------------------------------------------
T_steps        = 5000;
dt_env         = 0.02;
dt_tcm         = 1/600;
n_tcm_substeps = 10;

% Pong world parameters
env.width        = 1.0;
env.height       = 1.0;
env.paddle_x     = 0.1;
env.paddle_h     = 0.15;
env.paddle_speed = 0.6;
env.ball_speed   = 0.6;

% Actions
A_UP   = 1;
A_STAY = 2;
A_DOWN = 3;
actions = [A_UP, A_STAY, A_DOWN];

% Controller / learning parameters
ctrl.k_intercept        = 5.0;      % slope of target tanh for intercept error
ctrl.max_t_intercept    = 3.0;      % ignore absurdly far future intercepts
ctrl.approach_threshold = -1e-4;    % vx must be negative enough to count as approaching
ctrl.eta_W              = 2e-3;     % normalized LMS learning rate
ctrl.weight_decay       = 1e-5;     % mild stabilising decay on readout
ctrl.readout_clip       = 10.0;     % bound W_motor values
ctrl.target_gate_floor  = 0.05;     % weak learning when ball is not approaching
ctrl.plot_every         = 1;        % set to e.g. 2/5/10 if plotting is slow

% Action selection parameters
plan.horizon_steps      = 25;       % short-horizon rollout
plan.lambda_motor       = 0.25;     % cost for deviating from TCM motor prior
plan.lambda_effort      = 0.03;     % cost for moving rather than staying
plan.lambda_switch      = 0.04;     % cost for changing action
plan.lambda_hit_bonus   = 0.60;     % reward-like bonus for likely paddle contact
plan.discount           = 0.96;     % temporal discount in rollout

% -------------------------------------------------------------------------
% Initialise environment
% -------------------------------------------------------------------------
x_env = zeros(5,1);
x_env(1) = 0.8;
x_env(2) = 0.5;

theta0   = pi + (pi/3)*(rand - 0.5);   % broadly leftwards
x_env(3) = env.ball_speed * cos(theta0);
x_env(4) = env.ball_speed * sin(theta0);
x_env(5) = 0.5;

% -------------------------------------------------------------------------
% Initialise TCM
% -------------------------------------------------------------------------
[P_tcm, M_tcm, x_tcm0] = init_tcm_demo();
x_tcm = x_tcm0(:);

% Put TCM in a simple stable-ish initial state, as in original demo.
x_tcm = x_tcm * 0;
x_tcm(1:min(8,numel(x_tcm))) = -50;

Nx = numel(x_tcm);

% Sparse motor readout. Keep this aligned with useful states in your model.
idx_motor = [2 4 10 12];
idx_motor = idx_motor(idx_motor <= Nx);
if isempty(idx_motor)
    idx_motor = 1:min(8,Nx);
end

W_motor = zeros(1, Nx);
W_motor(idx_motor) = 1e-3 * randn(1, numel(idx_motor));

% -------------------------------------------------------------------------
% Logging
% -------------------------------------------------------------------------
traj.x_env         = zeros(5, T_steps);
traj.v_pred        = zeros(1, T_steps);
traj.v_target      = zeros(1, T_steps);
traj.dy_now        = zeros(1, T_steps);
traj.dy_intercept  = zeros(1, T_steps);
traj.y_intercept   = nan(1, T_steps);
traj.t_intercept   = nan(1, T_steps);
traj.approaching   = false(1, T_steps);
traj.action        = zeros(1, T_steps);
traj.hit           = false(1, T_steps);
traj.miss          = false(1, T_steps);
traj.W_norm        = zeros(1, T_steps);

% -------------------------------------------------------------------------
% Main simulation loop
% -------------------------------------------------------------------------
a_t = A_STAY;
prev_action = A_STAY;

for t = 1:T_steps

    % 1) Update true environment given current action.
    [x_env, event] = step_pong_env(x_env, a_t, env, dt_env);

    % 2) Compute self-supervised intercept target from game physics.
    intercept = compute_intercept_target(x_env, env, ctrl);

    % 3) Build sensory summary for TCM.
    %    Keep the first components close to the original demo, but include
    %    intercept-relevant information when M.ninputs supports it.
    s_t = build_sensory_summary(x_env, intercept);

    % 4) Drive TCM with sensory input.
    x_tcm = step_tcm(x_tcm, s_t, P_tcm, M_tcm, dt_tcm, n_tcm_substeps);
    x_tcm = real(x_tcm(:));

    % 5) TCM motor readout.
    v_pred = tanh(W_motor * x_tcm);

    % 6) Learn motor readout from intercept target, not teacher controller.
    %    Gate learning strongly when the ball is approaching the paddle, but
    %    allow a tiny amount of stabilising learning otherwise.
    v_target = intercept.v_target;
    gate     = intercept.target_gate;
    err      = v_target - v_pred;

    x_m = x_tcm(idx_motor);
    denom = 1e-6 + (x_m' * x_m);
    dW = ctrl.eta_W * gate * err * (x_m' / denom);

    W_motor(idx_motor) = (1 - ctrl.weight_decay) * W_motor(idx_motor) + dW;
    W_motor(idx_motor) = max(-ctrl.readout_clip, min(ctrl.readout_clip, W_motor(idx_motor)));

    % Recompute readout after the update for control.
    v_pred = tanh(W_motor * x_tcm);

    % 7) Select next action using short-horizon model-based rollout.
    a_t = select_action_pong_intercept_rollout(x_env, v_pred, prev_action, env, dt_env, actions, plan);
    prev_action = a_t;

    % 8) Live dashboard.
    if mod(t, ctrl.plot_every) == 0
        draw_pong_frame_intercept(x_env, env, v_pred, v_target, a_t, intercept, event, t*dt_env);
    end

    % 9) Log.
    traj.x_env(:,t)        = x_env;
    traj.v_pred(t)         = v_pred;
    traj.v_target(t)       = v_target;
    traj.dy_now(t)         = x_env(2) - x_env(5);
    traj.dy_intercept(t)   = intercept.dy_intercept;
    traj.y_intercept(t)    = intercept.y_intercept;
    traj.t_intercept(t)    = intercept.t_intercept;
    traj.approaching(t)    = intercept.approaching;
    traj.action(t)         = a_t;
    traj.hit(t)            = event.hit;
    traj.miss(t)           = event.miss;
    traj.W_norm(t)         = norm(W_motor(idx_motor));
end

plot_results_pong_tcm_intercept(traj, env, dt_env);

end

% =========================================================================
% Environment
% =========================================================================
function [x_env_next, event] = step_pong_env(x_env, a_t, env, dt)
% STEP_PONG_ENV  One-step update of classic Pong world.
%
% x_env = [ball_x; ball_y; ball_vx; ball_vy; paddle_y]

event.hit  = false;
event.miss = false;

bx  = x_env(1);
by  = x_env(2);
vx  = x_env(3);
vy  = x_env(4);
py  = x_env(5);

% Paddle update.
switch a_t
    case 1
        py = py + env.paddle_speed * dt;
    case 2
        % stay
    case 3
        py = py - env.paddle_speed * dt;
end
py = max(0, min(env.height, py));

% Ball update.
bx = bx + vx * dt;
by = by + vy * dt;

% Top/bottom reflections.
if by <= 0
    by = -by;
    vy = -vy;
elseif by >= env.height
    by = 2*env.height - by;
    vy = -vy;
end

% Right wall reflection.
if bx >= env.width
    bx = 2*env.width - bx;
    vx = -vx;
end

% Paddle contact / miss.
if bx <= env.paddle_x
    if by >= py - env.paddle_h && by <= py + env.paddle_h
        bx = 2*env.paddle_x - bx;
        vx = -vx;
        event.hit = true;
    else
        event.miss = true;
        bx = 0.9 * env.width;
        by = rand * env.height;
        theta = pi + (pi/2)*(rand - 0.5);
        vx = env.ball_speed * cos(theta);
        vy = env.ball_speed * sin(theta);
    end
end

x_env_next = [bx; by; vx; vy; py];
end

% =========================================================================
% Self-supervised intercept target
% =========================================================================
function intercept = compute_intercept_target(x_env, env, ctrl)
% COMPUTE_INTERCEPT_TARGET
% Predict where the ball will be vertically when it reaches the paddle x.
% This is the core self-supervised target.

bx = x_env(1);
by = x_env(2);
vx = x_env(3);
vy = x_env(4);
py = x_env(5);

approaching = vx < ctrl.approach_threshold && bx > env.paddle_x;

if approaching
    t_intercept = (env.paddle_x - bx) / vx;  % vx < 0, so this is positive
else
    t_intercept = inf;
end

valid = approaching && isfinite(t_intercept) && t_intercept >= 0 && t_intercept <= ctrl.max_t_intercept;

if valid
    y_raw = by + vy * t_intercept;
    y_intercept = reflect_y_periodic(y_raw, env.height);
    dy_intercept = y_intercept - py;
    target_gate = 1.0;
else
    % When the ball is moving away or too far in the future, use a weak
    % centre/position target rather than a hand-coded dy teacher.
    y_intercept = nan;
    t_intercept = nan;
    dy_intercept = 0;
    target_gate = ctrl.target_gate_floor;
end

v_target = tanh(ctrl.k_intercept * dy_intercept);

intercept.y_intercept   = y_intercept;
intercept.t_intercept   = t_intercept;
intercept.dy_intercept  = dy_intercept;
intercept.v_target      = v_target;
intercept.target_gate   = target_gate;
intercept.approaching   = approaching;
intercept.valid         = valid;
end

function y = reflect_y_periodic(y_raw, H)
% REFLECT_Y_PERIODIC  Reflect an unconstrained y through top/bottom walls.
% Equivalent to folding the real line onto [0,H] with mirror symmetry.

period = 2 * H;
y_mod = mod(y_raw, period);
if y_mod <= H
    y = y_mod;
else
    y = period - y_mod;
end
end

% =========================================================================
% Sensory input and TCM integration
% =========================================================================
function s_t = build_sensory_summary(x_env, intercept)
% BUILD_SENSORY_SUMMARY
% Sensory vector. If M.ninputs is small, step_tcm will truncate this vector.
%
% Components:
%   1 dy_now          current ball-paddle vertical error
%   2 ball_x          current ball x
%   3 ball_vy         current ball vertical velocity
%   4 ball_vx         current ball horizontal velocity
%   5 dy_intercept    predicted intercept error
%   6 approach flag   ball moving toward paddle

bx = x_env(1);
by = x_env(2);
vx = x_env(3);
vy = x_env(4);
py = x_env(5);

dy_now = by - py;

s_t = [dy_now; bx; vy; vx; intercept.dy_intercept; double(intercept.approaching)];
end

function x_tcm = step_tcm(x_tcm, s_t, P, M, dt, n_sub)
% STEP_TCM  Integrate TCM for a few micro-steps.

if isfield(M,'ninputs') && ~isempty(M.ninputs)
    n_inputs = M.ninputs;
else
    n_inputs = numel(s_t);
end

u = zeros(n_inputs,1);
k = min(numel(s_t), n_inputs);
u(1:k) = s_t(1:k);

x_tcm = x_tcm(:);

for i = 1:n_sub
    [dx,~,~] = atcm.tc_hilge2(x_tcm, u, P, M);
    dx = real(dx(:));
    x_tcm = x_tcm + dt * dx;

    % Safety check against numerical explosion.
    bad = ~isfinite(x_tcm);
    if any(bad)
        warning('TCM state became non-finite; resetting bad entries to zero.');
        x_tcm(bad) = 0;
    end
end
end

function [P, M, x] = init_tcm_demo()
% INIT_TCM_DEMO  Load TCM parameters from tcm_init.mat.

p = mfilename('fullpath');
[p,~,~] = fileparts(p);
load([p filesep 'tcm_init.mat']);

P = tcm_init.pE;
M = tcm_init.M;
x = M.x;

if ~isfield(M,'ninputs') || isempty(M.ninputs)
    M.ninputs = 3;
end
end

% =========================================================================
% Action selection with short-horizon rollout
% =========================================================================
function a_opt = select_action_pong_intercept_rollout(x_env, v_des, prev_action, env, dt_env, actions, plan)
% SELECT_ACTION_PONG_INTERCEPT_ROLLOUT
% Evaluate each immediate action by rolling the environment forward under
% that action for a short horizon. This makes action selection less myopic.

costs = zeros(numel(actions),1);

desired_dir = sign(v_des);

for i = 1:numel(actions)
    a0 = actions(i);
    x_roll = x_env;
    J = 0;
    disc = 1;
    hit_seen = false;

    for h = 1:plan.horizon_steps
        [x_roll, event] = step_pong_env_deterministic(x_roll, a0, env, dt_env);

        by = x_roll(2);
        py = x_roll(5);
        bx = x_roll(1);
        vx = x_roll(3);

        % Prioritise alignment when ball is approaching/near paddle.
        approach_weight = 0.25 + 0.75 * double(vx < 0);
        proximity_weight = 1 + 2 * max(0, 1 - abs(bx - env.paddle_x) / env.width);

        align_cost = approach_weight * proximity_weight * (by - py)^2;

        if event.hit
            hit_seen = true;
            align_cost = align_cost - plan.lambda_hit_bonus;
        elseif event.miss
            align_cost = align_cost + 2 * plan.lambda_hit_bonus;
        end

        J = J + disc * align_cost;
        disc = disc * plan.discount;
    end

    u_dir = action_to_dir(a0);
    motor_cost  = (u_dir - desired_dir)^2;
    effort_cost = abs(u_dir);
    switch_cost = double(a0 ~= prev_action);

    costs(i) = J + plan.lambda_motor  * motor_cost + ...
                   plan.lambda_effort * effort_cost + ...
                   plan.lambda_switch * switch_cost;

    if hit_seen
        costs(i) = costs(i) - 0.1 * plan.lambda_hit_bonus;
    end
end

[~, idx] = min(costs);
a_opt = actions(idx);
end

function [x_env_next, event] = step_pong_env_deterministic(x_env, a_t, env, dt)
% Deterministic rollout version: no random reset after a miss.
% Instead, apply a large miss event and reflect the ball, so action costs are
% comparable without injecting rollout noise.

event.hit  = false;
event.miss = false;

bx  = x_env(1);
by  = x_env(2);
vx  = x_env(3);
vy  = x_env(4);
py  = x_env(5);

switch a_t
    case 1
        py = py + env.paddle_speed * dt;
    case 2
    case 3
        py = py - env.paddle_speed * dt;
end
py = max(0, min(env.height, py));

bx = bx + vx * dt;
by = by + vy * dt;

if by <= 0
    by = -by;
    vy = -vy;
elseif by >= env.height
    by = 2*env.height - by;
    vy = -vy;
end

if bx >= env.width
    bx = 2*env.width - bx;
    vx = -vx;
end

if bx <= env.paddle_x
    if by >= py - env.paddle_h && by <= py + env.paddle_h
        bx = 2*env.paddle_x - bx;
        vx = -vx;
        event.hit = true;
    else
        % Mark miss, then reflect as a deterministic continuation.
        bx = 2*env.paddle_x - bx;
        vx = -vx;
        event.miss = true;
    end
end

x_env_next = [bx; by; vx; vy; py];
end

function u_dir = action_to_dir(a)
switch a
    case 1
        u_dir = +1;
    case 2
        u_dir = 0;
    case 3
        u_dir = -1;
    otherwise
        u_dir = 0;
end
end

% =========================================================================
% Plotting / animation
% =========================================================================
function draw_pong_frame_intercept(x_env, env, v_pred, v_target, a_t, intercept, event, t)
% DRAW_PONG_FRAME_INTERCEPT
% Live dashboard for self-supervised intercept controller.

persistent hFig hAx

bx = x_env(1);
by = x_env(2);
py = x_env(5);
px = env.paddle_x;

if isempty(hFig) || ~isvalid(hFig)
    hFig = figure('Color','w');
    set(hFig,'Position',[200 100 950 800]);

    hAx.game      = subplot(3,4,[1 2 5 6]);
    hAx.motor     = subplot(3,4,[3 4]);
    hAx.target    = subplot(3,4,[7 8]);
    hAx.action    = subplot(3,4,[9 10]);
    hAx.intercept = subplot(3,4,[11 12]);

    subplot(hAx.game); hold on;
    rectangle('Position',[0 0 env.width env.height], ...
              'EdgeColor',[0.6 0.6 0.6], 'LineWidth',1.5);
    hAx.ball = plot(bx, by, 'ko', 'MarkerFaceColor','k', 'MarkerSize',10);
    hAx.paddle = plot([px px], [py-env.paddle_h, py+env.paddle_h], ...
                      'r-', 'LineWidth',6);
    hAx.intercept_marker = plot(nan, nan, 'gx', 'MarkerSize',10, 'LineWidth',2);
    hAx.intercept_line = plot([px px], [0 env.height], 'g:', 'LineWidth',1.0);
    axis([0 env.width 0 env.height]); axis manual;
    set(gca,'YDir','normal'); xlabel('x'); ylabel('y');
    title('Pong-TCM with self-supervised intercept target');

    subplot(hAx.motor); hold on;
    hAx.motor_line = plot(nan, nan, 'b', 'LineWidth',1.3);
    ylabel('v_{TCM}'); title('TCM motor drive');

    subplot(hAx.target); hold on;
    hAx.target_line = plot(nan, nan, 'k', 'LineWidth',0.9);
    hAx.pred_line   = plot(nan, nan, 'r', 'LineWidth',1.2);
    ylabel('drive'); legend({'intercept target','TCM'}, 'Location','eastoutside');
    title('Self-supervised target vs TCM readout');

    subplot(hAx.action); hold on;
    hAx.action_line = plot(nan, nan, 'b', 'LineWidth',1.2);
    yticks([1 2 3]); yticklabels({'UP','STAY','DOWN'});
    ylim([0.5 3.5]); ylabel('action'); title('Actions');

    subplot(hAx.intercept); hold on;
    hAx.dy_line = plot(nan, nan, 'm', 'LineWidth',1.2);
    hAx.hit_line = plot(nan, nan, 'g.', 'MarkerSize',10);
    hAx.miss_line = plot(nan, nan, 'r.', 'MarkerSize',10);
    ylabel('dy_{intercept}'); title('Predicted intercept error and events');

    guidata(hFig, hAx);
end

hAx = guidata(hFig);

% Game panel.
set(hAx.ball, 'XData', bx, 'YData', by);
set(hAx.paddle, 'XData', [px px], 'YData', [py-env.paddle_h, py+env.paddle_h]);
if intercept.valid
    set(hAx.intercept_marker, 'XData', px, 'YData', intercept.y_intercept);
else
    set(hAx.intercept_marker, 'XData', nan, 'YData', nan);
end

% Motor drive.
x = get(hAx.motor_line,'XData'); y = get(hAx.motor_line,'YData');
set(hAx.motor_line,'XData',[x t],'YData',real([y v_pred]));

% Target vs readout.
x1 = get(hAx.target_line,'XData'); y1 = get(hAx.target_line,'YData');
x2 = get(hAx.pred_line,'XData');   y2 = get(hAx.pred_line,'YData');
set(hAx.target_line,'XData',[x1 t],'YData',real([y1 v_target]));
set(hAx.pred_line,  'XData',[x2 t],'YData',real([y2 v_pred]));

% Action trace.
xa = get(hAx.action_line,'XData'); ya = get(hAx.action_line,'YData');
set(hAx.action_line,'XData',[xa t],'YData',real([ya a_t]));

% Intercept dy and events.
xd = get(hAx.dy_line,'XData'); yd = get(hAx.dy_line,'YData');
set(hAx.dy_line,'XData',[xd t],'YData',real([yd intercept.dy_intercept]));

if event.hit
    xh = get(hAx.hit_line,'XData'); yh = get(hAx.hit_line,'YData');
    set(hAx.hit_line,'XData',[xh t],'YData',[yh intercept.dy_intercept]);
elseif event.miss
    xm = get(hAx.miss_line,'XData'); ym = get(hAx.miss_line,'YData');
    set(hAx.miss_line,'XData',[xm t],'YData',[ym intercept.dy_intercept]);
end

drawnow limitrate;
end

function plot_results_pong_tcm_intercept(traj, env, dt_env)

t = (0:size(traj.x_env,2)-1) * dt_env;

by = traj.x_env(2,:);
py = traj.x_env(5,:);

n_hits = sum(traj.hit);
n_miss = sum(traj.miss);
hit_rate = n_hits / max(1, n_hits + n_miss);

figure('Color','w');

subplot(5,1,1);
plot(t, by, 'LineWidth',1.3); hold on;
plot(t, py, 'LineWidth',1.3);
ylabel('y'); legend({'ball y','paddle y'}, 'Location','eastoutside');
title(sprintf('Pong-TCM intercept controller: hits=%d, misses=%d, hit rate=%.2f', n_hits, n_miss, hit_rate));

subplot(5,1,2);
plot(t, traj.y_intercept, 'LineWidth',1.2); hold on;
plot(t, py, 'LineWidth',1.2);
ylabel('y'); legend({'predicted intercept y','paddle y'}, 'Location','eastoutside');
title('Predicted intercept target');

subplot(5,1,3);
plot(t, traj.v_target, 'LineWidth',1.1); hold on;
plot(t, traj.v_pred, 'LineWidth',1.1);
ylabel('drive'); legend({'intercept target','TCM readout'}, 'Location','eastoutside');
title('Self-supervised target vs learned TCM motor readout');

subplot(5,1,4);
plot(t, traj.action, 'LineWidth',1.0);
yticks([1 2 3]); yticklabels({'UP','STAY','DOWN'});
ylabel('action'); title('Selected actions');

subplot(5,1,5);
plot(t, traj.W_norm, 'LineWidth',1.1);
xlabel('time (s)'); ylabel('||W||'); title('Motor readout norm');

end
