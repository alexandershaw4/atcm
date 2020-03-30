
% Run the integration
[y,w,s,~,~,pst,l,n,f,QD,Sp] = feval(DCM.M.IS,Qp,DCM.M,DCM.xU);

% collect firing patterns
times     = Sp{1}(:,1);
CellFired = Sp{1}(:,2);

nc = unique(CellFired); % num unique cells (that fired!)

% sort firing into a matrix - cell-by-time
for i  = 1:length(nc)
    it = find(CellFired==nc(i)); 
    CF(i,:)         = pst*0;
    CF(i,times(it)) = 1;
end

% get the membrane potentials over time
mV = squeeze(s{1}(1,:,1,:));

for i = 1:8
    subplot(8,1,i);
    
    [r(i),p(i)] = corr( CF(i,:)' , mV(i,:)' );
    
    scatter( CF(i,:) , mV(i,:) ); lsline;
    title(sprintf('R^2 = %d , p = %d',r(i).^2,p(i)));
end