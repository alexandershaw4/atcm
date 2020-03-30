function [Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, GrowthRates, POD_Mode_Energies]=dmd(Big_X, r, dt)
%
%This function is a wrapper for performing DMD on reduced order models
%of fluid dynamical systems. The wrapper is built such that the common
%DMD outputs are easily available, so only the birds-eye theoretical
%grasp is needed for operation.

%====INPUTS===
%Big_X: The ND matrix to analyze the DMD modes. First dimension is
%time. All other dimensions are arbitrarily chosen and will be reshaped
%back when the eigenvectors are calculated.

%r: Number of reduced order POD modes to keep. Usually much less than
%the number of snapshots. If in doubt, process the data set with this function and
%plot stem(POD_Mode_Energies). The first few modes should have more
%energy than the remaining modes.

%dt: Time delay between snapshots. Gives the Nyquist limit of the data
%set and allows the output ModeFrequencies to be meaningful physical or
%dimensionless frequencies. If unused, make dt=1.

%====OUTPUTS===
%Eigenvalues: DMD mode complex eigenvalues (r by 1). Think that for each time
%step the corresponding eigenvector is multiplied by
%exp(Eigenvalue*dt). To understand more, google "z-transform".

%Eigenvectors: The complex mode shapes themselves (r by [N-1]D).
%Enables one to visualize the important dynamics of the system. Output
%is already in the same dimensions as the input, where the first (time) dimension now is
%replaced by mode number.

%ModeAmplitudes: The complex amplitudes of the mode shapes. Determines
%the dominant modes. i.e., ModeAmplitudes=max(abs(ModeAmplitudes)) is
%the dominant mode.

%ModeFrequencies: The real, two-sided frequencies related to each mode.
%DMD being a linear decomposition will have each mode be related to a
%single frequency.

%GrowthRates: The growth/decay rates related to each mode. Units [1/s].
%If GrowthRates<0, mode decays. If GrowthRates>0, mode grows. In
%realistic systems, modes GrowthRates<<0 will decay fast whereas modes
%with GrowthRates~0 probably correspond to a limit-cycle of some sort.

%POD_Mode_Energies: The mode energies of all POD modes. If the user
%doesn't know the appropriate value of the input "r", plotting
%stem(POD_Mode_Energies) or
%stem(cumsum(POD_Mode_Energies)/sum(POD_Mode_Energies)) can give
%insight on how many modes to keep. Good guesses would be to retain at
%least 70% of the total energy of the system.


dims=size(Big_X);
newDims=dims;
newDims(1)=r;
%Removes mean. Note: Not removing the mean biases the modes as the
%data points centroid is shifted. If one wants to capture only the
%oscillations around the mean, the mean MUST be removed.
Big_X=Big_X-repmat(mean(Big_X,1),[dims(1) ones(1,length(dims)-1)]);

%Reshapes Big_X
Big_X=(reshape(Big_X,dims(1),prod(dims(2:end)))).';


%Split Big_X into two snapshot sets
X=Big_X(:,1:end-1);
Y=Big_X(:,2:end);

%SVD on X
[U, S, V]=svd(X,0);

%Before reducing rank returns the mode energies for further analysis of
%the ROM validity
POD_Mode_Energies=diag(S).^2;

%Reduce rank
U=U(:,1:r);
V=V(:,1:r);
S=S(1:r,1:r);

%Gets A_tilde
A_tilde=U'*Y*V/S;

%(For debugging), we can compare if A_tilde=A, for r=max(r):
% A=Y*pinv(X);

%Compute A_tilde eigenvalues and eigenvectors
[eVecs, Eigenvalues] = eig(A_tilde);

%Gets the DMD eigenvectors back
Eigenvectors=Y*V*inv(S)*eVecs;
Eigenvalues=diag(Eigenvalues);

%Gets the mode amplitudes
ModeAmplitudes=Eigenvectors\X(:,1);

%Gets the frequencies associated with the modes
fNY=1/(2*dt);
ModeFrequencies=(angle(Eigenvalues)/pi)*fNY;

%Gets the growth rates
GrowthRates=log(abs(Eigenvalues))/dt;

%Reshapes the Eigenvectors back to original Big_X dims
Eigenvectors=Eigenvectors.';
Eigenvectors=reshape(Eigenvectors,newDims);

end