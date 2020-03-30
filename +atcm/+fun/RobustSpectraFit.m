function [mnewspectra,Freqs,uncorrected,modelm,newspectra,FitPar ]=RobustSpectraFit(Freqs,singlespectra,PostSmooth)
f0=findthenearest(1,Freqs);
f1=findthenearest(95,Freqs);

if Freqs(f0)==0,
    f0=f0+1;
end

Freqs=log(Freqs(f0:f1));

fc1=atcm.fun.findthenearest(1.5,Freqs);
fc2=atcm.fun.findthenearest(4.25,Freqs);
FitFreqs=[Freqs(1:fc1) Freqs(fc2:end)];

N=size(singlespectra,1);

for j=1:N,
    m=singlespectra(j,:);
    lm=log(m(f0:f1));
    Fitlm=[lm(1:fc1) lm(fc2:end)];

    b=robustfit(FitFreqs,Fitlm);
    FitPar(j,:)=b;
    modelm(j,:)=(b(1)+b(2)*Freqs);
    
    uncorrected(j,:)=lm;
    
    lm=lm-modelm(j,:);
    lm=atcm.fun.moving_average(lm,PostSmooth);
    newspectra(j,:)=(lm);
end
mnewspectra=mean(newspectra,1);

end

