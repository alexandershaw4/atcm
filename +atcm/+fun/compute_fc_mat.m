function Z = compute_fc_mat(D,band,fs)

Dcat = reshape(D,[size(D,1) size(D,2)*size(D,3)]);

% filter -
if nargin > 1 && length(band) == 2 && nargin == 3
    for i = 1:size(Dcat,1)
        Dcat(i,:) = bandpassfilter(Dcat(i,:),fs,band);
    end
end

for i = 1:size(Dcat,1)
    Dcat(i,:) = abs(hilbert(Dcat(i,:)));
end

for i = 1:size(Dcat,1)
    for j = 1:size(Dcat,1)
        [r(i,j),p(i,j)] = corr( D(i,:)',D(j,:)' );
    end
end

Z=kzscorenan(r);

end

function Z=kzscorenan(Z)
s=nanstd(Z(:));
m=nanmean(Z(:));
Z=((Z-m)/s);
end