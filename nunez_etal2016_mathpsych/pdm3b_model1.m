function pdm3b_model1(inputfile,rmtrials,varargin)
%PDM3B_MODEL1 - Runs a new JAGS model without EEG inputs
%
%load jagsins.mat to see structure
%
%Usage: pdm3b_model1('jagsins.mat',rmtrials);
%
%
%Inputs:
%  inputfile: name of file that contains .eeg structure with subject
%                level EEG fields (i.e. jagsins.mat)
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
%  8/26/15        Michael Nunez                 Original code
%  3/04/16        Michael Nunez               Added thinning parameter that
%             exists in the original fits for the paper, but forgotten here

%% Initial

data = load(inputfile);

%Organize date and time into a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});

modelname = timestr;
cleanup = true;
nsamples = 5e3;
nburnin = 2e3;
thin = 10;
nchains =6;
verbosity =1;
parallelit = 1;
maxcores = 3;
modules = {'wiener' 'dic'};

%% JAGS code for the diffusion model

model = {
    'model {'
      '# No effect on b (bias between responses)'
    '# b <- .5'
    '# Boundary separation kept constant'
    '# a <- 1'
    ''
    '# Effect on T_er (preprocessing /nondecision time)'
    '# Varies by subject'
    'tsd ~ dgamma(5, 20)' %x=linspace(0, 1, 100); plot(x, gampdf(x, 5, .05))
    'ttau <- pow(tsd, -2)'
    'for (c2 in 1:3) {'  %noise
    '    tmu[c2] ~ dnorm(.3, 1)T(0,3)' %std = 1
    '    for (sub in 1:nsubs) {'  %subject
    '       t[c2,sub] ~ dnorm(tmu[c2], ttau)'
    '    }'
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
    '# Effect on v (drift rate)'
    '# Varies by subject'
    'vsd ~ dgamma(5, 5)' % x=linspace(0, 4, 100); plot(x, gampdf(x, 5, .2))
    'vtau <- pow(vsd, -2)'
    'for (c2 in 1:3) {'  %noise
    '    vmu[c2] ~ dnorm(1.5, 1/16)T(-9,9)' %std = 4
    '    for (sub in 1:nsubs) {'  %subject
    '       v[c2,sub] ~ dnorm(vmu[c2], vtau)'
    '    }'
    '}'
    ''
    '# Likelihood'
    'for (i in 1:n)'
    '{'
    ''
    '    y[i] ~ dwiener(1/s[noise[i],subject[i]], t[noise[i],subject[i]], 0.5, v[noise[i],subject[i]]/s[noise[i],subject[i]])'
    '}'
    '}'
    };

%% Code for Trinity

params = {'tsd' 'tmu' 't' 'ssd' 'smu' 's' 'vsd' 'vmu' 'v'};

samprt = data.rt;
rt = samprt/1024; %Reaction time in samples to seconds (not milliseconds)
ntrials = length(data.correct);

%Remove no answer trials and subjects that were taken out of this model
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
    'v', randn(3,R.nsubs));


%Setup data split index for cross-validation
rmindex = logical(rmtrials(~tremove));

%Training data
S.nsubs = R.nsubs;
S.subject = R.subject(~rmindex);
S.n = R.n - sum(rmindex);
S.y = R.y(~rmindex);
S.jitter = R.jitter(~rmindex);
S.noise = R.noise(~rmindex);

%Test data
T.nsubs = R.nsubs;
T.subject = R.subject(rmindex);
T.n = R.n - sum(~rmindex);
T.y = R.y(rmindex);
T.jitter = R.jitter(rmindex);
T.noise = R.noise(rmindex);

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
    'cutoffs','info','params','S','T');

