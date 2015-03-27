function COV = spatialcov(x,XY,D,normalize)
%
% Compute sample autocovariance of field variable measured at a set of
% spatial sites as a function of lag distance.
%
% Input:
%  x            = variable of interest
%  XY           = (2 x N) array of site coordinates (in a local Cartesian system)
%  D            = (d x 1) vector of lag distance bins
%  normalize    = 1: normalize by sample variance of x (i.e., autocorrelation 
%                 coefficient function), 0: otherwise.
%
% Output:
%   COV         = (d x 1) vector of sample autocovariance as a function of lag distance
%
% Andreas Mavrommatis 2015

Nsta = length(XY);

mu = mean(x);

COV = nan(length(D)-1,1);
for k = 1:length(D)-1
    Nk = 0; C = 0;
    for i = 1:Nsta
        
        % compute distances from each reference station to all test stations
        ref_xy = XY(:,i);                % reference station coordinates
        d = nan(Nsta,1);
        for j = 1:Nsta
            test_xy = XY(:,j);           % test station coordinates
            d(j) = norm(ref_xy - test_xy);  % distance
        end
        
        % find and append number of stations at or within within each distance
        keep = find(d < D(k+1) & d >= D(k));
        Nk = Nk + length(keep);
        
        % Compute and append (unnormalized) sample covariance between
        % each reference station and all its test stations
        C = C + sum( (x(i) - mu)*(x(keep) - mu) );
    end
    
    % Sample covariance for each lag distance
    COV(k) = C/Nk;
    
end

% Normalize by sample variance?
if normalize
    COV = COV/var(x);
end
