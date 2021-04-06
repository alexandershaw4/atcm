function GenMeanDataset(list,name)

Data  = [];
for i = 1:length(list)
    clear D;
    D    = spm_eeg_load(list{i});
    Data = cat(3,Data,D(:,:,:));
end

name = [name '_MeanDataset'];
MeanDataset = clone(D,name,[size(Data)]);
MeanDataset(:,:,:) = Data;
MeanDataset.save;