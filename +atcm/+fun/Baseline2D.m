 function Output = Baseline2D(Input, Start, Finish, Type)
 %Usage Output = Baseline2D(Input, Start, Finish, Type)
 %Data must be in Channels * Samples format or Frequencies x Samples
 %Type is a dummy variable; if present a relative change Baseline will be calculated
 Dims = size(Input);
 NSamps = Dims(2);
 Baseline = mean(Input(:,Start:Finish),2);
 A = repmat(Baseline, [1, NSamps]);
 if nargin == 3 %This is a subtracted Baseline 
     
     Output = Input - A;
 elseif nargin == 4 %This is the relative CHange Baseline
     Output = 100 * (Input -  A) ./ A;  
 end
     