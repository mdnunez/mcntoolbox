function [cutoff,dmatewmaplot]=ewmav2(data,L,l,s)
% EWMAV2  Determine RT cut-off with the EWMA method
%   [CUTOFF, DMATEWMAPLOT] = EWMAV2(DATA,LIMITS,LAMBDA,SIGMA), where DATA
%   is a regular N-by-3 data matrix, and LIMITS, LAMBDA, and SIGMA are
%   parameters for the Exponentially Weighted Moving Average algorithm as
%   described in the DMAT Manual.
%   CUTOFF is the cut-off value below which RT values are censored.
%   DMATEWMAPLOT is a structure that can be used as input for PLOTEWMA.
%
%   See also PLOTEWMA, OUTLIERTREATMENT.
%
%   Author: Joachim Vandekerckhove (joachim.vandekerckhove@psy.kuleuven.be)
%   Part of the DMA Toolbox. Please read the End User License Agreement,
%   contained in 'dmateula.txt' or by invoking the DMATLICENSE command. 
%   See also http://ppw.kuleuven.be/okp/dmatoolbox.

% Note: If you like one-liners, sort the data into XX. Then,
%   cutoff=T(find(diff(filter(l,[1 (l-1)],XX,(1-l)/2)'<.5+L*s*...
%           sqrt(l/(2-l)*(1-(1-l).^(2*(1:size(XX,1))))))==-1,1,'first'));
% jv

%% Check input
if ~isvaliddataset(data,3) && ~isvaliddataset(data,2)
    error('DMAT:ewmav2:invalidDataSet',...
        'Data set is not valid Nx3 data matrix.')
end
if ~isscalar(l)
    error('DMAT:ewmav2:incorrectInput','''Lambda'' should be a scalar.')
end
if ~isscalar(L)
    error('DMAT:ewmav2:incorrectInput','''Limits'' should be a scalar.')
end
if ~isscalar(s)
    error('DMAT:ewmav2:incorrectInput','''Sigma'' should be a scalar.')
end

%% Sort data
[T,I]=sort(data(:,end));
X=data(I,end-1);

%% Define constants
cutoff=0;

%% Apply EWMA filter
ucl=.5+L*s*sqrt(l/(2-l)*(1-(1-l).^(2*(1:size(X,1)))));
lcl=.5-L*s*sqrt(l/(2-l)*(1-(1-l).^(2*(1:size(X,1)))));
z=filter(l,[1 (l-1)],X,(1-l)/2)';
vv=find(diff(z<ucl)==-1,1,'first');

if ~isempty(vv)
    cutoff=T(vv)+eps;
else
    vv = 1;
end

%% Finish output
T=T';
dmatewmaplot=struct('RTs',T,...
    'UCL',ucl,'LCL',lcl,...
    'ControlState',s,...
    'BlueX',[cutoff T(T>cutoff)],...
    'BlueY',[ucl(vv) z(T>cutoff)],...
    'RedX',[T(T<cutoff) cutoff],...
    'RedY',[z(T<cutoff) ucl(vv)],...
    'Lambda',l,'Limits',L);
