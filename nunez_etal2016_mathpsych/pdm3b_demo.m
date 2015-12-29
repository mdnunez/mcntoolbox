%% README

%THIS SCRIPT WILL NOT WORK WITHOUT SPECIFIC PACKAGES AND PROGRAMS
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
%  12/29/15        Michael Nunez                 Original code

%% Training/Test Split

%2/3 - 1/3 Training/Test split
oneszeros1 = [zeros(1,360) ones(1,180)];
oneszeros2 = [zeros(1,320) ones(1,160)];

rng('default');
rng(11);
splitdata1 = [];
for j=1:17
    splitdata1 = [splitdata1 oneszeros1(randperm(length(oneszeros1)))];
end

splitdata2 = [];
for j=1:3
    splitdata2 = [splitdata2 oneszeros2(randperm(length(oneszeros2)))];
end

whichsubs1 = 1:17;
%Note that randperm's use changed after a specific MATLAB version
randsubs1 = whichsubs1(randperm(17,4));

for j = randsubs1
    splitdata1([1:540] + 540*(j-1)) = ones(1,540);
end

%% Code

pdm3b_model1('jagsins.mat','nsamples',5e3,'nburnin',2e3,...
    'thin',10,'nchains',6,'maxcores',3,'rmtrials',splitdata1);

pdm3b_model2('jagsins.mat',{'p200trialerp_c1n' 'p200trialerplat_c1n' 'n200trialerp_c1r' 'n200trialerplat_c1r'},...
    'nsamples',5e3,'nburnin',2e3,...
    'thin',10,'nchains',6,'maxcores',3,'rmtrials',splitdata1);

pdm3b_model3('jagsins.mat',{'p200trialerp_c1n' 'p200trialerplat_c1n' 'n200trialerp_c1r' 'n200trialerplat_c1r'},...
    'nsamples',5e3,'nburnin',2e3,...
    'thin',10,'nchains',6,'maxcores',3,'rmtrials',splitdata1);