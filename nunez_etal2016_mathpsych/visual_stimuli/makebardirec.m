function barimage = makebardirec(imagesize,barsize,imagetype,jitter,preload, graphit)
%makebardirec - Builds a bar patch oriented in two ways
%
%   This function creates a matrix of bars with mean rotation directed
%   either north-east or north-west on average
%
%Useage: 
%  >> barimage = makebardirec(imagesize,density,imagetype,jitter,preload,graphit);
%
%Useage Example: 
%  >> barimage = makebardirec(600,12,1,0,0,1);
%
%To test try:
%  >> a = randi(2,1); b = rand(1)*90; makebardirec(600,12,a,b,0,1);
%
%Inputs:
%   imagesize = number of pixels for each side of the image (pixels = imagesize*imagesize) (e.g., 600)
%
%   density = density parameter (use of 6 in pdmexp2)
%
%   imagetype (optional): default is random
%   1 - lines directed north-east / south-west on average
%   2 - lines directed north-west / south-east on average
%
%   jitter (optional): default is 0, which is a figure with 'exact' rotations
%   Acceptable values on the range [0 120], corresponding to possible
%   degrees of random rotation
%
%   preload (optional): default is 0, 1 calls preloaded/(45/135)_(jitter).mat
%   file to load the bar image maps
%
%   graphit (optional): Graphs the image in grayscale if = 1

%% README

%For use with pdmexp3_public

%% Possible Citations
% Nunez, M. D., Srinivasan, R., & Vandekerckhove, J. (2015). Individual differences in attention influence perceptual decision making. Frontiers in psychology, 8.
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
%  
%% Record of revisions:
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  5/20/13      Michael Nunez                    Original code
%  1/1/16       Michael Nunez       Defaults better reflect stimuli used in experiments              

%% Arguments

switch nargin
    case 0
    imagesize = 600;
    barsize = 6;
    imagetype = randsample(1:2,1);
    jitter = 60; 
    preload = 0;
    graphit = 1;
    case 1
    barsize = 6;
    imagetype = randsample(1:2,1);
    jitter = 60;
    preload = 0;
    graphit = 1;
    case 2
    imagetype = randsample(1:2,1);
    jitter = 60;
    preload = 0;
    graphit = 1;
    case 3
    jitter = 60;
    preload = 0;
    graphit = 1;
    case 4
    preload = 0;
    graphit = 1;
    case 5
    graphit = 0;
end

if isempty(imagetype)
    imagetype = randsample(1:2,1);
end;

%% Load Preloaded Bars 
if preload
    if imagetype == 1
        anglefilename = '45';
    else
        anglefilename = '135';
    end
    if exist(['preloaded/',anglefilename,'_',num2str(jitter),'.mat']) == 2
        load(['preloaded/',anglefilename,'_',num2str(jitter),'.mat']);
        eval(['loadedbars = m',anglefilename,'j',num2str(jitter),';']);
    else
        error('These specific bar image maps do not exist in the subdirectory ''preloaded''.');
    end
end

%% Useful pre-code calculations
halfimage = round(imagesize/2);
halfbar = floor(barsize/2);
halfbarwidth = round(barsize/12);
halfbarlength = floor(barsize/2);


%Density calculation
density = 3*barsize; %Density calculation

%Jitter to radians
jitterrad = jitter*(pi/180); 

%% Circular Display
circlech = ones(imagesize);
xc = 1:imagesize;
yc = 1:imagesize;
z = nan(imagesize);
for i = 1:imagesize
    for j = 1:imagesize
        z(i,j) = sqrt((xc(i)-halfimage)^2 + (yc(j)-halfimage)^2);
        if z(i,j) > round(5*imagesize/12);
            circlech(i,j) = 0;
        end
    end
end

%% Build a unrotated bar
%patch = zeros(barsize);
patch = zeros(12);
%Contrust bar with length of barsize and width
%2*round(barsize/12)
%patch((halfbar - halfbarlength+1):(halfbar + halfbarlength),(halfbar - halfbarwidth+1):(halfbar + halfbarwidth))= 1;
patch(1:12,6:7) = 1;

%% Build an oriented bar patch
barimage = zeros(imagesize);


[xg,yg]=meshgrid(-halfbar:halfbar,-halfbar:halfbar);


%Draw the image with bar patches
for i= 1:(round(imagesize/density)-1) %No longer pixel by pixel (as is in the egg paradigm)
    for j = 1:(round(imagesize/density)-1)
        %Put some random placement into the bars
        randplacex = randi(round(1.5*barsize),1) - round(.75*barsize);
        randplacey = randi(round(1.5*barsize),1) - round(.75*barsize);

        %Find the correct indicies
        wherey = i*density + randplacey;
        wherex = j*density + randplacex;
        if all(all(circlech((wherey-halfbar):(wherey+halfbar),(wherex-halfbar):(wherex+halfbar))))
            
            %Here we insert our rand jitter to the figure angles
            randjitter = (2*rand(1)-1)*jitterrad;

            %Radians to angles
            randjitter = randjitter*(180/pi);
            if imagetype == 1
                figangle(i,j) = -45 + randjitter;
            elseif imagetype == 2
                figangle(i,j) = 45 + randjitter;
            end
            
            if ~preload
                %Draw the figure-bar
                rotpatch = imrotate(patch, figangle(i,j));
            else
                %Draw from prerotated images
                randindex = randi(length(loadedbars),1);
                rotpatch = loadedbars{randindex};
            end
            
            %Find size of rotpatch
            [row column] = size(rotpatch);
            
            %Place the figure-bar into the bar patch
            barimage((wherey-floor(row/2)+1):(wherey+ceil(row/2)),(wherex-floor(column/2)+1):(wherex+ceil(row/2))) = rotpatch;
        end
    end
end

barimage(barimage == 1) = .5;

%% Graph the image

if graphit == 1
%figure;
imagesc(barimage,[0 1]);
colormap(gray);
screensize = get(0,'ScreenSize');
set(gcf,'Position', [1 1 screensize(3) screensize(4)]);
axis('square');
end