function evalERPplot(DCM,P,G)

pst = DCM.xY.pst;
Y   = DCM.xY.y;
y   = atcm.fun.evalERP(DCM,P,G);

for i = 1:length(y)
    subplot(1,length(y),i);
    plot(pst,Y{i},':',pst,y{i});
end