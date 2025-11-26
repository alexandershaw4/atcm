function W = make_visual_gamma_mask(t, f)
W = ones(numel(f), numel(t));
% Emphasise post-onset 0.1–1.2 s
W(:, t<0.1 | t>1.2) = 0.25;
% Emphasise 30–90 Hz band (gamma)
W(f<30 | f>90, :) = W(f<30 | f>90, :) * 0.3;
% De-emphasise <10 Hz (big slow ERPs)
W(f<10,:) = W(f<10,:)*0.2;
end
