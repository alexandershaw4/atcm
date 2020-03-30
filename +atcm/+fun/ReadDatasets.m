function [Datasets] = ReadDatasets(listfile)
%
%Reads in all the Datasets specified by listfile and checks that they exist
%
fid = fopen(listfile, 'r');
if fid < 0, error('cannot open listfile'); end

i = 0;
%fprintf('Datasets to process...\n\n');
while feof(fid) == 0
  filler =  fgetl(fid);
  deblank(filler);
  if length(filler) == 0 
  
  else    
    i = i + 1;
    Datasets(i) = cellstr(filler);
    %fprintf('     %s\n',filler);
    Dataset = cell2mat(Datasets(i));
%     if exist(Dataset) ==0    %Check that each of the datasets exists before continuing
%        error(['Could not find the dataset " ' Dataset '"']);
%     end
    if Datasets{i}(end) == '/' %Remove tailing slashes accidentally included
       Datasets{i}(end) = [];
    end   
  end  
end 
fclose(fid);
