%pdm3b_trialsortplot - Script to create single-trial evoked responses plot, Figure 5 of Nunez et al., 2017 Journal of Mathematical Psychology 
%
% Copyright (C) 2016 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  9/28/15       Michael Nunez                   Original code
%  12/15/15      Michael Nunez                Save high quality image
%   2/26/16      Michael Nunez                   Color bar

%% Sort single-trial EEG

load('sublist.mat');
prepdata = load('svdins3b_5.mat');
load('svdinputs8.mat');

noisetrials = nan(768,540*17);
resptrials = nan(1024,540*17);
noisetrialsrt = nan(768,540*17);
resptrialsrt =nan(1024,540*17);
for s=1:17
    thesetrials = (1:540) + 540*(s-1);
    
    %Sort by max magnitude
    [~,tempindx] = sort(prepdata.p200trialerp_c1n(prepdata.subject == s));
    noisetrials(:,thesetrials) = squeeze(trialerp{1,s}(:,1,tempindx));
    [~,tempindx2] = sort(prepdata.p200trialerp_c1r(prepdata.subject == s));
    resptrials(:,thesetrials) = squeeze(trialerp{2,s}(:,1,tempindx2));
    
%     %Sort by min magnitude
%     [~,tempindxrth] = sort(prepdata.rt(prepdata.subject == s & prepdata.noise == .6),2,'descend');
%     [~,tempindxrtm] = sort(prepdata.rt(prepdata.subject == s & prepdata.noise == .45),2,'descend');
%     [~,tempindxrtl] = sort(prepdata.rt(prepdata.subject == s & prepdata.noise == .3),2,'descend');
%     tempindxrt = [tempindxrth tempindxrtm tempindxrtl];
%     noisetrialsrt(:,thesetrials) = squeeze(trialerp{1,s}(:,1,tempindxrt));
%     resptrialsrt(:,thesetrials) = squeeze(trialerp{2,s}(:,1,tempindxrt));
    
    %Organize by RT tertiles magnitude
    noises = [.6 .45 .3];
    for n=1:3
    temprt = prepdata.rt(prepdata.subject == s & prepdata.noise == noises(n));
%     tert = prctile(temprt,[33 67]);
%     noisetrialsrt(:,(1:180)+180*(n-1) + 540*(s-1)) = ...
%         [noisetrials(:,(temprt > tert(2))) noisetrials(:,(temprt <= tert(2)) & (temprt > tert(1))) ...
%         noisetrials(:,(temprt <= tert(1))) noisetrials(:,isnan(temprt))];
%     resptrialsrt(:,(1:180)+180*(n-1) + 540*(s-1)) = ...
%         [resptrials(:,(temprt > tert(2))) resptrials(:,(temprt <= tert(2)) & (temprt > tert(1))) ...
%         resptrials(:,(temprt <= tert(1))) resptrials(:,isnan(temprt))];
    rtmedian = nanmedian(temprt);
    noisetrialsrt(:,(1:180)+180*(n-1) + 540*(s-1)) = ...
        [noisetrials(:,(temprt > rtmedian)) ...
        noisetrials(:,(temprt <= rtmedian)) noisetrials(:,isnan(temprt))];
    resptrialsrt(:,(1:180)+180*(n-1) + 540*(s-1)) = ...
        [resptrials(:,(temprt > rtmedian)) ...
        resptrials(:,(temprt <= rtmedian)) resptrials(:,isnan(temprt))];
    end
end


%%Find p200 and n200 locations and magnitudes

n200 = round(150*1.024):round(275*1.024);
p200 = round(150*1.024):round(275*1.024);

noisetrials = flipud(noisetrials');
resptrials = flipud(resptrials');
noisetrialsrt = flipud(noisetrialsrt');
resptrialsrt = flipud(resptrialsrt');

noisedelaysrt = zeros(size(noisetrialsrt));
respdelaysrt = zeros(size(resptrialsrt));
[~,nind] = max(noisetrialsrt(:,p200),[],2);
[~,rind] = min(-resptrialsrt(:,n200),[],2);
for n=1:540
    noisedelaysrt(n,nind(n)) = 1;
    respdelaysrt(n,rind(n)) = 1;
end

%% Plot single-trial EEG for one subject

windowsize = round(300*1.024);

%%Paper figure
whichsub = 12;
thisindx = (1:540) + (whichsub-1)*540;

xms = [0:50:300];
xsamps = round(xms*1.024);

fontsize = 20;

f1 = figure('units','normalized','outerposition',[0 .5 1 .5]);
% f1 = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1);
%  subplot(2,2,1);
%sub1 = imagesc(noisetrials(thisindx,1:windowsize),[-100 100]);
sub1 = imagesc(noisetrials(thisindx,1:windowsize),[-90 90]);
set(gca,'XLim',[0 windowsize],'XTick',xsamps,'XTickLabel',xms,'YTickLabel',[],'Fontsize',16);
ylabel('Cue Intervals','Fontsize',fontsize);
xlabel('Delay (ms) following noise onset','Fontsize',fontsize);
line(round([150 150]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
line(round([275 275]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
%line(round([p200loc(whichsub) p200loc(whichsub)]*1.024),[1 540],'Color',[1 1 1],'LineWidth',2,'Linestyle','--');
%line(round([p200loc2(whichsub) p200loc2(whichsub)]*1.024),[1 540],'Color',[0 1 0],'LineWidth',2,'Linestyle','--');

subplot(1,2,2);
%  subplot(2,2,2);
%sub2 = imagesc(-resptrials(thisindx,1:windowsize),[-100 100]);
sub2 = imagesc(-resptrials(thisindx,1:windowsize),[-90 90]);
set(gca,'XLim',[0 windowsize],'XTick',xsamps,'XTickLabel',xms,'YTickLabel',[],'Fontsize',16);
ylabel('Response Intervals','Fontsize',fontsize);
xlabel('Delay (ms) following signal onset','Fontsize',fontsize);
line(round([150 150]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
line(round([275 275]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
%line(round([n200loc(whichsub) n200loc(whichsub)]*1.024),[1 540],'Color',[1 1 1],'LineWidth',2,'Linestyle','--');
%line(round([n200loc2(whichsub) n200loc2(whichsub)]*1.024),[1 540],'Color',[0 1 0],'LineWidth',2,'Linestyle','--');
sub2ax = gca;
originalsize = get(sub2ax,'Position');
cb = colorbar('EastOutside','Fontsize',fontsize,'YTick',[-90 -45 0 45 90]);
yl = ylabel(cb,'\muV','Fontsize',fontsize);
set(yl,'Rotation',270);
set(sub2ax,'Position',originalsize);


% subplot(2,2,3);
% %sub3 = imagesc(noisetrialsrt(thisindx,1:windowsize),[-100 100]);
% sub3 = imagesc(noisetrialsrt(thisindx,1:windowsize),[-90 90]);
% %sub3 = imagesc(noisedelaysrt(thisindx,1:windowsize));
% set(gca,'XLim',[0 windowsize],'XTick',xsamps,'XTickLabel',xms,'YTickLabel',[],'Fontsize',16);
% ylabel('Noise Intervals','Fontsize',fontsize);
% xlabel('Delay (ms) following noise onset','Fontsize',fontsize);
% line(round([150 150]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
% line(round([275 275]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
% % line(round([p200loc(whichsub) p200loc(whichsub)]*1.024),[1 540],'Color',[1 1 1],'LineWidth',2,'Linestyle','--');
% % line(round([p200loc2(whichsub) p200loc2(whichsub)]*1.024),[1 540],'Color',[0 1 0],'LineWidth',2,'Linestyle','--');
% templine = downsample(nind(thisindx),10)+153;
% %line(templine,10:10:540,'Color',[0 0 0],'LineWidth',1.5);
% line([0 windowsize], [181 181],'Color','k','LineWidth',2);
% line([0 windowsize], [361 361],'Color','k','LineWidth',2);
% 
% subplot(2,2,4);
% %sub4 = imagesc(-resptrialsrt(thisindx,1:windowsize),[-100 100]);
% sub4 = imagesc(-resptrialsrt(thisindx,1:windowsize),[-90 90]);
% set(gca,'XLim',[0 windowsize],'XTick',xsamps,'XTickLabel',xms,'YTickLabel',[],'Fontsize',16);
% ylabel('Response Intervals','Fontsize',fontsize);
% xlabel('Delay (ms) following signal onset','Fontsize',fontsize);
% line(round([150 150]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
% line(round([275 275]*1.024),[1 540],'Color','k','LineWidth',2,'Linestyle','--');
% % line(round([n200loc(whichsub) n200loc(whichsub)]*1.024),[1 540],'Color',[1 1 1],'LineWidth',2,'Linestyle','--');
% % line(round([n200loc2(whichsub) n200loc2(whichsub)]*1.024),[1 540],'Color',[0 1 0],'LineWidth',2,'Linestyle','--');
% templine = downsample(rind(thisindx),10) + 153;
% %line(templine,10:10:540,'Color',[0 0 0],'LineWidth',1.5);
% line([0 windowsize], [181 181],'Color','k','LineWidth',2);
% line([0 windowsize], [361 361],'Color','k','LineWidth',2);

%% Save figure

export_fig(f1,'ERPs_single_trials','-opengl','-png','-r200');