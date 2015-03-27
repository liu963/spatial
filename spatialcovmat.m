function Sigma = spatialcovmat(sampvar,staxy,model,dodisp)

% Populate empirical covariance matrix of spatial data given a model of distance-
% based autocovariance.
% 
% Input: 
%   sampvar = sample variance of field variable (scalar)
%   staxy   = (2 x N) array of site coordinates
%   model   = fit object after modeling autocovariance as a function of distance
%   dodisp  = 1: for test plot, 0: otherwise
%
% Output:
%   Sigma   = (N x N) covariance matrix.
%
% Andreas Mavrommatis, 2015


if (nargin < 4), dodisp = false; end

n = length(staxy);
   
% Form distance matrix of stations
D = dist(staxy);

% Evaluate covariances based on given sample variance, model, and distances
F = sampvar*feval(model,D);

% Shape into square matrix
Sigma = reshape(F,n,n);


% Testing
if dodisp
    
    % plot matrix
    figure, imagesc(Sigma), colorbar
    title('\Sigma')
    
    % plot a random row 
    i = randi(length(Sigma));  
    c = Sigma(i,:);
    d = nan(n,1);
    ref_xy = staxy(:,i);
    for j = 1:n
        test_xy = staxy(:,j);           % test station coordinates
        d(j) = norm(ref_xy - test_xy);  % distance from reference station to test station
    end
    cmax = max(c);
    
    ds = linspace(min(d),max(d));
    fv = feval(model,ds);
    c_pred = sampvar*fv;
    
    
    figure, hold on
    plot(d(c~=cmax),c(c~=cmax), '.')
    plot(ds,c_pred,'r')
    xlabel('distance from station (km)')
    ylabel('covariance')
    title(['Row ', num2str(i)])

    

end

end
