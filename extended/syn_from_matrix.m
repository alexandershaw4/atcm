function Syn = syn_from_matrix(M, pop, mix)
% mix: optional struct of receptor proportion splits
%   mix.exc = struct('AMPA',0.6,'NMDA',0.3,'KA',0.1);
%   mix.inh = struct('GABAA',0.8,'GABAB',0.2);

if nargin < 3
  mix.exc = struct('AMPA',0.6,'NMDA',0.3,'KA',0.1);
  mix.inh = struct('GABAA',0.8,'GABAB',0.2);
end

exc_pre = {'SS','SP','DP','TP','RC'};
inh_pre = {'SI','DI','RT'};
cortical = {'SS','SP','SI','DP','DI','TP'};
thalamic = {'RT','RC'};

rows = {};
for to = 1:numel(pop)
  for from = 1:numel(pop)
    w = M(to,from);
    if w <= 0, continue; end
    pre  = pop{from};  post = pop{to};

    if ismember(pre, inh_pre)
      receptors = fieldnames(mix.inh)';  % {'GABAA','GABAB'}
      pre_mods  = {}; if ismember(pre,{'SI','DI'}), pre_mods = {'MOR'}; end
    else
      receptors = fieldnames(mix.exc)';  % {'AMPA','NMDA','KA'}
      pre_mods  = {}; if strcmp(pre,'SP'), pre_mods = {'CB1'}; end
    end

    if     ismember(pre,cortical) && ismember(post,cortical), delay = 'Intra';
    elseif ismember(pre,cortical) && ismember(post,thalamic), delay = 'Cx2Th';
    else                                                      delay = 'Th2Cx';
    end

    post_mods = {}; if any(strcmp(post,{'SP','DP'})), post_mods = {'5HT2A'}; end
    stp = double(strcmp(pre,'SP'));

    for r = 1:numel(receptors)
        rec = receptors{r};

        % choose the right mixture group based on presyn class
        if ismember(pre, inh_pre)
            group = mix.inh;     % e.g. fields: GABAA, GABAB
        else
            group = mix.exc;     % e.g. fields: AMPA, NMDA, (KA) ...
        end

        % robust lookup (skip if receptor not present in chosen mix)
        if ~isfield(group, rec)
            % either skip this receptor or set frac=0
            % continue
            frac = 0;  % <- alternative: keep row with zero weight
        else
            frac = group.(rec);
        end

        base_wr = w * frac;
        param   = sprintf('w_%s_%s_%s', pre, post, rec);
        rows(end+1,:) = {pre,post,rec,delay,stp,pre_mods,post_mods,param,base_wr}; %#ok<AGROW>
    end  
  end
end

Syn = cell2table(rows, 'VariableNames', ...
 {'pre','post','receptor','delay_family','stp_id','pre_mods','post_mods','param','base_w'});
end
