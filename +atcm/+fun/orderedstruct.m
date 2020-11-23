function new = orderedstruct(obj1,obj2)

f = fieldnames(obj1);

for i = 1:length(f)
    
    try
        new.(f{i}) = obj2.(f{i});
    catch
        new.(f{i}) = spm_unvec( spm_vec(obj1.(f{i}))*0, obj1.(f{i}) );
    end
    
end