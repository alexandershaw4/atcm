function g = gate_value(channel_name, V, P)
switch channel_name
  case 'NMDA'
    s = 0.2; a = 0.062; sc = 1; % scale params; add to P if needed
    g = 1 ./ (1 + s * exp(-a * sc * V));
  otherwise
    g = 1;
end
end
