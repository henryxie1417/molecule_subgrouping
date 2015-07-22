function  [Msdout,SDarray3D,SDarray2D,SDarray1D] = MSDalldim(track,dim)
%-------------------------------------------------------------------------
%
% EXAMPLE:
%       [Msdout,SDarray3D,SDarray2D,SDarray1D] = MSDalldim(track,dim);
%
% PURPOSE:
%       calculate the mean-square displacement and its standard deviation/
%       variance for each particle. the latter are calculated according to
%       Quian et al, BJ 60 (1991) 910
%
% INPUTS:
%       tracksin = track matrix [framenumber, X, Y, Z]
%       dim      = dimension of msd analysis
%
% OUTPUTS:
%       Msdout    = tracknum, framenum, MSD, dMSD, nSD 
%       SDarrayXD = SD arrays for all dimensions
%
% CALLS:
%       nothing
%
% HISTORY: 
%       adapted from MSD_Ver7.m from Thomas Schmidt, Gerhard Blab, and 
%       Laurent Holtzer (Leiden University), last changes as of 
%       20070118 (v4.00) by LH
%
%       20080626 TOB    first version, changed from original to
%                       only do one track at a time 
%           
% TODO: 
%       nothing
%
%--------------------------------------------------------------------------

%% input checks
if nargin<1, help MSDalldim, return, end
if nargin<2, dim=3; end

%% preallocate some variables
Msdout    = [];
SDarray3D = {};
SDarray2D = {};
SDarray1D = {};
if isempty(track), return, end


%% MSD from 1 to Nlag-1
Nlag = size(track,1); % number of points in this track
L1 = tril(repmat((2:Nlag)',1,Nlag-1));
L2 = fliplr(triu(repmat((1:Nlag-1)',1,Nlag-1)));
L1 = L1(any(L1(:),2));
L2 = L2(any(L2(:),2));
L  = track(L1,1)-track(L2,1);
H1D = (track(L2,2:4)-track(L1,2:4)).^2;
H2D = sum(H1D(:,1:2),2);
H3D = sum(H1D,2);
 
%% calculate SD arrays if necessary
Mlag = track(end,1)-track(1,1); %length of track in images
if size(SDarray3D,2)<Mlag
    [SDarray3D{size(SDarray3D,2)+1:Mlag}] = deal([]); %create empty arrays
    [SDarray2D{size(SDarray2D,2)+1:Mlag}] = deal([]); %of right size
    [SDarray1D{size(SDarray1D,2)+1:Mlag}] = deal([]);
end

if nargout>1
    SDarray3Dtemp = accumarray(L, H3D, [], @(x) {x})';
    SDarray2Dtemp = accumarray(L, H2D, [], @(x) {x})';
    L3 = repmat(L,1,3);
    L3(:,2) = L3(:,1)+max(L(:,1));
    L3(:,3) = L3(:,2)+max(L(:,1));
    SDarray1Dx = accumarray(L3(:),H1D(:),[], @(x) {x})';
    for i=1:length(SDarray1Dx)/3
        SDarray1Dtemp{i} = [SDarray1Dx{i},SDarray1Dx{i+length(SDarray1Dx)/3}, ...
            SDarray1Dx{i+2*length(SDarray1Dx)/3}];
    end
    for i=1:max(L)
        SDarray3D{i} = [SDarray3D{i}; SDarray3Dtemp{i}];
        SDarray2D{i} = [SDarray2D{i}; SDarray2Dtemp{i}];
        SDarray1D{i} = [SDarray1D{i}; SDarray1Dtemp{i}];
    end
end

%% calculate the msd dependent on dim
if dim==1,     H = H1D(:,1);
elseif dim==2, H = H2D;
else,          H = H3D;
end
Lag = unique(L);
Msd = accumarray(L,H,[],@mean);
Msd = Msd(Lag);
nSD = histc(L,Lag)';

%% calculate variance in msd according to Quian et al, BJ 60 (1991) 910
dMsd = zeros(size(Msd));
dMsd(nSD==1) = Msd(nSD==1);

%[il,K,sqrt((4*il^2*K+2*K+il-il^3)/(6*il*K^2))]
ind = nSD>=Lag';
K   = nSD(ind);
il  = Lag(ind)';
dMsd(ind) = Msd(ind) .* (sqrt((4*il.^2.*K+2*K+il-il.^3)./(6*il.*K.^2)))';

%[il,K,sqrt(1+(K^3-4*il*K^2+4*il-K)/(6*il^2*K))]
ind = nSD<Lag';
K  = nSD(ind);
il = Lag(ind)';
dMsd(ind) = Msd(ind) .* (sqrt(1+(K.^3-4*il.*K.^2+4*il-K)./(6*il.^2.*K)))';

%% MSD
Msdout = [Lag,Msd,dMsd,nSD'];

