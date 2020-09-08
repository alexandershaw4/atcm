function names = convertnames_htoconn(hn)
% convert shorthand intrinsic para name to full - 
% e.g. h6 --> 'ss -> tp'
%      hn64 --> 'nmda rl -> rl'
%
%  names = convertnames_htoconn(hn)
%
% AS

C = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};

n = 0;
for i = 1:length(C)
    for j = 1:length(C)
        n = n + 1;
        list{n} = sprintf('%s -> %s',C{i},C{j});
    end
end

list = reshape(list,[8 8]);

id = strrep(hn,'H','');
id = strrep(id,'n','');

nmda = find(contains(hn,'n'));

for i = 1:length(id)
    num(i) = str2num(id{i});
end

names = list(num);

for i = 1:length(nmda)
    names{nmda(i)} = ['nmda ' names{nmda(i)}];
end