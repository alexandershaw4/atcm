function [w_eff, STP] = stp_step(w_base, m_pre, STP, dt)
% m_pre: presyn firing per edge (vector)
% simple rate-based surrogate
du = (STP.u*0 + STP.u) .* (1 - STP.u)/STP.F + STP.u .* m_pre; %#ok<NASGU>
dR = (1 - STP.R)/STP.D - STP.u .* STP.R .* m_pre; %#ok<NASGU>
w_eff = w_base .* STP.u .* STP.R; % placeholder
end
