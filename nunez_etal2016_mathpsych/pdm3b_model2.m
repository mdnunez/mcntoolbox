function pdm3b_model2(inputfile,rmtrials,eegfields,varargin)
%PDM3B_MODEL2 - Runs a new JAGS Model Type 2 with EEG inputs
%
%load jagsins.mat to see structure
%
%Usage: pdm3b_model2('jagsins.mat',rmtrials,eegfields);
%
%
%Inputs:
%  inputfile: name of file that contains .eeg structure with subject
%                level EEG fields (i.e. jagsins.mat)
%  eegfields: EEG fields of (eeginputfile) to include in model
%
%% README

%THIS PROGRAM WILL NOT WORK WITHOUT SPECIFIC PACKAGES AND PROGRAMS
%
%JAGS - Just Another Gibbs Sampler
%http://sourceforge.net/projects/mcmc-jags/
%
%jags-wiener - Wiener distribution functions for JAGS
%http://sourceforge.net/projects/jags-wiener/
%
%DMAT - Diffusion Model Analysis Toolbox
%https://ppw.kuleuven.be/okp/software/dmat/
%
%Trinity
%https://github.com/joachimvandekerckhove/trinity
%
%Others?

%% Citation
% Nunez M.D., Vandekerchove, J., Srinivasan, R. (2016) How attention influences perceptual decision making: Single-trial EEG correlates of drift-diffusion model parameters. Journal of Mathematical Psychology.

%% Copyright 2015 Michael D. Nunez

%This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  8/31/15        Michael Nunez                 Original code
%  3/04/16        Michael Nunez               Added thinning parameter that
%             exists in the original fits for the paper, but forgotten here

%% Initial

data = load(inputfile);

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

%EEG effect on model parameters?
for t=1:length(eegfields)
    defaulttsv{t} = [1 1 1];
end

tsv = defaulttsv;
modelname = timestr;
cleanup = true;
nsamples = 5e3;
nburnin = 2e3;
nchains =6;
thin = 10;
verbosity =1;
parallelit = 1;
maxcores = 3;
modules = {'wiener' 'dic'};


%% JAGS code for the diffusion model

%Set up string cell array with proper indexing
%Set up cell structure for creating mu equations
tarray = [];
sarray = [];
varray = [];
for f=1:numel(eegfields)
    neweegfields{f} = strrep(eegfields{f},'_',''); %JAGS does not like the underscore in a data variable name
    covariates{f} = [neweegfields{f} '[i]']; %Subscripts for the .jags code
    if tsv{f}(1)
        tarray = cat(2,tarray,[covariates(f) num2cell(f)]'); 
    end
    if tsv{f}(2)
        sarray = cat(2,sarray,[covariates(f) num2cell(f)]');
    end
    if tsv{f}(3)
        varray = cat(2,varray,[covariates(f) num2cell(f)]');
    end
end

model = {
    'model {'
    '# No effect on b (bias between responses)'
    '# b <- .5'
    '# Boundary separation kept constant'
    '# a <- 1'
    ''
    '# 1 unit increase of eegfield{f} is associated with this additive effect on t'
    sprintf('for (f in 1:%i) {',size(tarray,2))
    'tbetasd[f] ~ dgamma(5, .2)' % x=linspace(0, 100, 100); plot(x, gampdf(x, 5, 5))
    'tbetatau[f] <- pow(tbetasd[f], -2)'
    '   for (c2 in 1:3) {' %noise
    '      tbetamu[c2,f] ~ dnorm(0,.0001)' %std = 100
    '      for (sub in 1:nsubs) {'  
    '         tbeta[c2,sub,f] ~ dnorm(tbetamu[c2,f],tbetatau[f])' %std = 100, These should be drawn from a EEG parameter specific distribution with some std parameter
    '      }'
    '    }'
    '}'
    ''
    '# Effect on Ter (non-decision time)'
    '# Varies by subject'
    'tsd ~ dgamma(5, .2)' % x=linspace(0, 100, 100); plot(x, gampdf(x, 5, 5))
    'ttau <- pow(tsd, -2)'
    'for (c2 in 1:3) {'  %noise
    '   tmu[c2] ~ dnorm(0,.0001)' %std = 100
    '   for (sub in 1:nsubs) {'  
    '      talpha[c2,sub] ~ dnorm(tmu[c2], ttau)'
    '   }'
    '}'
    ''
    '# Effect on s (diffusion coefficient)'
    '# Varies by subject'
    'ssd ~ dgamma(5, 20)' % x=linspace(0, 1, 100); plot(x, gampdf(x, 5, .05))
    'stau <- pow(ssd, -2)'
    'for (c2 in 1:3) {'  %noise
    '    smu[c2] ~ dnorm(.6, 1/4)T(0,4)' %std = 2
    '    for (sub in 1:nsubs) {'  %subject
    '       s[c2,sub] ~ dnorm(smu[c2], stau)'
    '    }'
    '}'
    ''
    '# 1 unit increase of eegfield{f} is associated with this additive effect on v'
    sprintf('for (f in 1:%i) {',size(tarray,2))
    'vbetasd[f] ~ dgamma(5, .2)' % x=linspace(0, 100, 100); plot(x, gampdf(x, 5, 5))
    'vbetatau[f] <- pow(vbetasd[f], -2)'
    '   for (c2 in 1:3) {'  %noise
    '      vbetamu[c2,f] ~ dnorm(0,.0001)' %std = 100
    '      for (sub in 1:nsubs) {'  %subject
    '         vbeta[c2,sub,f] ~ dnorm(vbetamu[c2,f],vbetatau[f])' %std = 100, These should be drawn from a EEG parameter specific distribution with some std parameter
    '      }'
    '    }'
    '}'
    ''
    '# Effect on v (diffusion process between trials)'
    '# Varies by condition and subject'
    'vsd ~ dgamma(5, .2)' % x=linspace(0, 100, 100); plot(x, gampdf(x, 5, 5))
    'vtau <- pow(vsd, -2)'
    'for (c2 in 1:3) {'  %noise
    '   vmu[c2] ~ dnorm(0,.0001)' %std = 100
    '   for (sub in 1:nsubs) {'
    '       valpha[c2,sub] ~ dnorm(vmu[c2], vtau)'
    '   }'
    '}'
    ''
    '# Likelihood'
    'for (i in 1:n)'
    '{'
    ''
    sprintf('v[i] <- valpha[noise[i],subject[i]]%s',sprintf(' + %s*vbeta[noise[i],subject[i],%i]',sarray{:}))
    ''
    ''
    sprintf('t[i] <- talpha[noise[i],subject[i]]%s',sprintf(' + %s*tbeta[noise[i],subject[i],%i]',sarray{:}))
    ''
    '    y[i] ~ dwiener(1/s[noise[i],subject[i]], t[i], 0.5, v[i]/s[noise[i],subject[i]])'
    '}'
    '}'
    };

%% Code for Trinity

params = {'tsd' 'tmu' 'talpha' 'tbeta' 'tbetamu' 'tbetasd' 'ssd' 'smu' 's' ...
    'vsd' 'vmu' 'valpha' 'vbeta' 'vbetamu' 'vbetasd'};

samprt = data.rt;
rt = samprt/1024; %Reaction time in samples to seconds (not milliseconds)
ntrials = length(data.correct);

%Remove no answer trials and subjects that weretaken out of this model
tremove = isnan(data.correct) | (data.goodtrials == 0); %Remove this last flag if imputing bad trials

%Remove RTs less than cutoff given by ewmav2
for j=1:length(data.subname)
    [cutoff] = ewmav2([data.correct(data.subject == j & ~isnan(rt))' rt(data.subject == j & ~isnan(rt))'],2,.01,.5);
    tremove = tremove | (rt < cutoff & data.subject == j);
    cutoffs(j) = cutoff;
end
nremove = sum(tremove);

correct = data.correct;
correct(correct == 0) = -1;
y = correct(~tremove).*rt(~tremove);

R.subject = data.subject(~tremove);
R.subject = R.subject(:);

tempc = data.noise;
tempc = tempc(~tremove);
[~, ~, noise] = unique(tempc);

tempj = data.jitter;
tempj = tempj(~tremove);
[~, ~, jitter] = unique(tempj);

R.n = ntrials - nremove;
R.nsubs = length(data.subname); %Keep this fixed for the total number of subjects even if removing subjects from modeling
R.y = y(:);
R.jitter = jitter(:);
R.noise = noise(:);

initstruct = @()struct(...
    'valpha', randn(3,R.nsubs),'talpha', .1 + rand(3,R.nsubs)*.2);


for f=1:numel(eegfields)
    R.(neweegfields{f}) = data.(eegfields{f})(~tremove)';
end


%Setup data split index for cross-validation
rmindex = logical(rmtrials(~tremove));

%Training data
S.nsubs = R.nsubs;
S.subject = R.subject(~rmindex);
S.n = R.n - sum(rmindex);
S.y = R.y(~rmindex);
S.jitter = R.jitter(~rmindex);
S.noise = R.noise(~rmindex);
for f=1:numel(eegfields)
    S.(neweegfields{f}) = R.(neweegfields{f})(~rmindex);
end

%Test data
T.nsubs = R.nsubs;
T.subject = R.subject(rmindex);
T.n = R.n - sum(~rmindex);
T.y = R.y(rmindex);
T.jitter = R.jitter(rmindex);
T.noise = R.noise(rmindex);
for f=1:numel(eegfields)
    T.(neweegfields{f}) = R.(neweegfields{f})(rmindex);
end

%% Run JAGS

fprintf('Building JAGS model %s and saving output...',modelname);

tic
[stats, chains, diagnostics, info] = callbayes('jags', ...
    'model', model, ...
    'data', S, ...
    'cleanup', cleanup, ...
    'nsamples', nsamples, ...
    'nburnin', nburnin, ...
    'nchains', nchains, ...
    'thin',thin,...
    'verbosity', verbosity, ...
    'workingdir','wdir', ...
    'monitorparams', params, ...
    'parallel',parallelit, ...
    'maxcores',maxcores, ...
    'modules',modules, ...
    'init', initstruct, ...
    varargin{:}); 

info.comptime = toc/60;
fprintf('JAGS took %f minutes!\n', info.comptime)

save(sprintf('jagsmodel%s.mat',modelname),'stats', 'chains', 'diagnostics',...
    'cutoffs','info','params','S','T','tsv');

