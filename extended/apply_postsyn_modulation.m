function gain = apply_postsyn_modulation(channel_name, Reg, Comp, P, ns, np)
mods = fieldnames(Comp.postMap);
gain = ones(ns,np);
for m = 1:numel(mods)
  mname = mods{m};
  if ~isfield(P, mname), continue; end
  sens = exp(P.(mname));
  % If any edges for this channel are tagged with this postsyn mod, scale gain
  hits = strcmp(Reg.Syn.receptor, channel_name) & Comp.postMap.(mname);
  if any(hits)
    gain = gain .* (1 + 0.2*sens); % example postsyn potentiation
  end
end
end
