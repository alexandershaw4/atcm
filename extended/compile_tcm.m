function Comp = compile_tcm(Reg, ns)
%COMPILE_TCM  Build sparse operators and index maps from registries.
%  ns: number of sources (areas)
% Returns Comp with fields:
%  .idx    — maps names->indices (pops, channels, mods, delays)
%  .A      — struct of receptor->sparse operator (per-source block diag)
%  .delay  — sparse delay matrix by families
%  .stp    — STP state containers per edge family (optional)
%  .postMap, .preMap — modulation selector matrices

np = numel(Reg.Pops);    
% Index maps
idx.pop = containers.Map({Reg.Pops.name}, num2cell(1:np));
idx.chan = containers.Map({Reg.Channels.name}, num2cell(1:numel(Reg.Channels)));
idx.mod  = containers.Map({Reg.Mods.name},    num2cell(1:numel(Reg.Mods)));
idx.delay= containers.Map({Reg.Delays.name},  num2cell(1:numel(Reg.Delays)));

% Receptor→channel map (same names here)
rec2chan = idx.chan;

% Build per-receptor adjacency (within each source; extend to cross-source if needed)
receptors = unique(Reg.Syn.receptor);
for r = 1:numel(receptors)
  thisR = receptors{r};
  A = sparse(np*ns, np*ns);
  rows = find(strcmp(Reg.Syn.receptor,thisR));
  for k = rows'
    i_pre  = idx.pop(Reg.Syn.pre{k});
    i_post = idx.pop(Reg.Syn.post{k});
    for s = 1:ns
      lin_pre  = (s-1)*np + i_pre;
      lin_post = (s-1)*np + i_post;
      A(lin_post, lin_pre) = 1; % weight applied later via params
    end
  end
  Comp.A.(thisR) = A;
end

% Delay matrix (family blocks)
N = np*ns; Delay = sparse(N,N);
for k = 1:height(Reg.Syn)
  i_pre  = idx.pop(Reg.Syn.pre{k});
  i_post = idx.pop(Reg.Syn.post{k});
  d_ms   = Reg.Delays(idx.delay(Reg.Syn.delay_family{k})).value_ms;
  for s = 1:ns
    lin_pre  = (s-1)*np + i_pre;
    lin_post = (s-1)*np + i_post;
    Delay(lin_post, lin_pre) = d_ms/1000; % seconds
  end
end
Comp.delay = Delay;

% Modulation selector maps (edge-target lists)
Comp.preMap  = build_mod_map(Reg, idx, 'pre_mods');
Comp.postMap = build_mod_map(Reg, idx, 'post_mods');

% Channel indices (where conductances live in x-state): supply externally or assume order
Comp.idx = idx;
Comp.rec2chan = rec2chan;
end

function ModMap = build_mod_map(Reg, idx, field)
np = numel(Reg.Pops); mods = {Reg.Mods.name};
Nedges = height(Reg.Syn);
ModMap = struct();
for m = 1:numel(mods)
  name = mods{m}; sel = false(Nedges,1);
  for k = 1:Nedges
    sel(k) = any(strcmp(Reg.Syn.(field){k}, name));
  end
  ModMap.(name) = sel; % boolean mask over edges of Syn
end
end
