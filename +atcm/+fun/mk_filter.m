function [B, A] = mk_filter(Fs,hpHz, lpHz, FilterOrder)

 
  Wn=[hpHz/Fs lpHz/Fs]*2;
  if Wn(1) <= 0,
      Wn(1)=2/Fs;
  end
  if Wn(2) >= 1,
      Wn(2)=1-(2/Fs);
  end
  
  [B,A] = butter(FilterOrder,Wn);
  return