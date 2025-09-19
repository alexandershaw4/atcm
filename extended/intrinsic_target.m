function T = intrinsic_target(ch, V, gNMDA, P)
switch ch
  case 'M'    % Kv7: more open with depolarization
    Vh = -35; k = 8;     T = 1 ./ (1 + exp(-(V - Vh)/k));
  case 'HCN'  % Ih: more open with hyperpolarization
    Vh = -80; k = 6;     T = 1 ./ (1 + exp((V - Vh)/k));
  case 'SK'   % Ca-activated K+: proxy Ca from NMDA
    alpha = 1.0; Ca = alpha * gNMDA; T = Ca ./ (1 + Ca);  % bounded 0..1
  case 'NaP'  % persistent Na+: depolarization-activated window
    Vh = -55; k = 7;     T = 1 ./ (1 + exp(-(V - Vh)/k));
  otherwise
    T = 0;
end
end
