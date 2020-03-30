function [Px,BestCost] = abcAS(fun,Px,pC)
%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA114
% Project Title: Implementation of Artificial Bee Colony in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%


% Problem Definition
CostFunction = fun;        % Cost Function

nVar      = length(Px);     % Number of Decision Variables
VarSize   = [1 nVar];       % Decision Variables Matrix Size
VarMin    = -1;             % Decision Variables Lower Bound
VarMax    = +1;             % Decision Variables Upper Bound

MaxIt     = 100;            % Maximum Number of Iterations
nPop      = 50;%100;            % Population Size (Colony Size)
nOnlooker = nPop;           % Number of Onlooker Bees
L = round(0.6*nVar*nPop);   % Abandonment Limit Parameter (Trial Limit)
a = 1;                      % Acceleration Coefficient Upper Bound

pC = spm_vec(pC)';

% Initialization
empty_bee.Position = [];        % Empty Bee Structure
empty_bee.Cost = [];

pop = repmat(empty_bee,nPop,1); % Initialize Population Array

BestSol.Cost = Inf;        % Initialize Best Solution Ever Found
BestCost = zeros(MaxIt,1); % Array to Hold Best Cost Values
C = zeros(nPop,1);         % Abandonment Counter

fprintf('Initial Bee Population Search...\n');
for i = 1:nPop % Create Initial Population
    phi1 = unifrnd(VarMin,VarMax,VarSize);
    pop(i).Position = Px'+(pC.*phi1); % TA: pC ??
    pop(i).Cost = CostFunction(pop(i).Position);
    if pop(i).Cost <= BestSol.Cost
        BestSol = pop(i);
    end
end

% ABC Main Loop
fprintf('Initialising main optim loop...\n');
for it = 1:MaxIt
    fprintf('Iteration %d / %d\n',it,MaxIt)
    
    % Recruited Bees
    fprintf('(1.) Recruiting\n');
    for i = 1:nPop
        K = [1:i-1 i+1:nPop]; % Choose k randomly, not equal to i
        k = K(randi([1 numel(K)]));
        phi = a*pC.*unifrnd(VarMin,VarMax,VarSize); % Define Acceleration Coeff.
        newbee.Position = pop(i).Position+phi.*(pop(i).Position-pop(k).Position); % New Bee Position
        newbee.Cost = CostFunction(newbee.Position); % Evaluation
        if newbee.Cost <= pop(i).Cost % Comparision
            pop(i) = newbee;
        else, C(i) = C(i)+1;
        end
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F = zeros(nPop,1);
    
    % alex: de-nan
    for i = 1:nPop
        if isnan( pop(i).Cost )
            pop(i).Cost = 0 ;
        end
    end
    
    MeanCost = mean([pop.Cost]);
    for i = 1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P = F/sum(F);
    
    % Onlooker Bees
    fprintf('(2.) Onlookers\n');
    for m = 1:nOnlooker
        i = atcm.optim.RouletteWheelSelection(P); % Select Source Site
        K = [1:i-1 i+1:nPop]; % Choose k randomly, not equal to i
        k = K(randi([1 numel(K)]));
        phi = a*pC.*unifrnd(-1,+1,VarSize); % Define Acceleration Coeff.
        newbee.Position = pop(i).Position+phi.*(pop(i).Position-pop(k).Position);% New Bee Position
        newbee.Cost = CostFunction(newbee.Position); % Evaluation
        if newbee.Cost <= pop(i).Cost % Comparision
            pop(i) = newbee;
        else, C(i) = C(i)+1;
        end
    end
    
    % Scout Bees
    fprintf('(3.) Scouts\n');
    for i = 1:nPop
        if C(i) >= L
            pop(i).Position = unifrnd(VarMin,VarMax,VarSize);
            pop(i).Cost = CostFunction(pop(i).Position);
            C(i) = 0;
        end
    end
    
    % Update Best Solution Ever Found
    for i = 1:nPop
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end
    
    BestCost(it) = BestSol.Cost; % Store Best Cost Ever Found
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]); % Display Iteration Information
    Px = BestSol.Position;
    
end


end

