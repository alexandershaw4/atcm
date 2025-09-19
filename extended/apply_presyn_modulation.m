function Drive = apply_presyn_modulation(Drive, Reg, Comp, P)
% Scale per-receptor Drive according to presynaptic modulators that target edges
mods = fieldnames(Comp.preMap);
for m = 1:numel(mods)
  mname = mods{m};
  if ~isfield(P, mname), continue; end
  sens = exp(P.(mname)); % sensitivity
  mask = Comp.preMap.(mname); % logical over edges
  % Simple implementation: scale the whole Drive for receptors that have edges under this mod
  for r = fieldnames(Drive)'
    rec = r{1}; %#ok<FXSET>
    if any(strcmp(Reg.Syn.receptor(mask), rec))
      Drive.(rec) = Drive.(rec) .* (1 ./ (1 + sens)); % suppress release
    end
  end
end
end

