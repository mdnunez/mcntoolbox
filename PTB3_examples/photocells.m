function photocells
%% README

%THIS PROGRAM WILL NOT WORK WITHOUT SPECIFIC CHANGES FOR YOUR SOFTWARE AND HARDWARE

%This program serves as examples of steady-state visual evoked
%potential (SSVEP) stimuli. See the two papers for a description of example
%experimental paradigms 
%
%In order to replicate this experimental stimulus, you will need to run 
%Psychtoolbox 3 on Linux with a good video card, and then test 
%the frequency of the "photocell" squares in each corner of the screen to 
%ensure the intended frequencies of the flickering stimuli and thus possible SSVEP
%responses (SSVEPs could also depend upon subject behavior, brain state,
%contrast, luminance, etc. and must be tested with preliminary studies)

%% Possible Citations
% Nunez, M. D., Srinivasan, R., & Vandekerckhove, J. (2015). Individual differences in attention influence perceptual decision making. Frontiers in psychology, 8.
% Nunez, M. D., Vandekerckhove, J., & Srinivasan, R. (2017) How attention influences perceptual decision making: Single-trial EEG correlates of drift-diffusion model parameters. Journal of Mathematical Psychology, 76, 117-130.

%% Copyright 2016 Michael D. Nunez

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
%  03/09/2016     Michael Nunez                 Original code

%% Initial Code

PsychJavaTrouble; 

AssertOpenGL; %Issue warning if PTB3 with non-openGL used

if ~IsLinux
    error('This program was written to run on Ubuntu Linux.');
end

%% Experimenter Prompt

%Inputs Prompt and Output Setup
%Experimenter Prompt
Screenres = get(0,'Screensize');

prompt1={'Window Pointer:',...
'Screen Length (x-axis):','Screen Width (y-axis):','Refresh Rate (fps):',...
'Flicker1 Freq (Hz):','Flicker2 Freq (Hz):','Number of Trials:'};
def1={'0',num2str(Screenres(3)),num2str(Screenres(4)),'120','40','30','10'};
studytitle='Photocell Test';


lineNo=1;
answer=inputdlg(prompt1,studytitle,lineNo,def1);

%Window Pointer / 'Home Screen'.  0 - the primary monitor; 1 - the secondary monitor.
whichScreen = str2num(answer{1});
%Screen resolution on the x-axis
xres = str2num(answer{2});

%Screen resolution on the y-axis
yres = str2num(answer{3});

%This should be the same as the Refresh Rate shown in the Display
%Properties on the computer.  Always check before running the experiment to
%match flicker frequency.
%This code is currently set up to only handle multiples of 60 fps.
refrate = str2num(answer{4});
realrefrate = Screen(0,'FrameRate');
if refrate ~= Screen(0,'FrameRate')
    error(['The real screen refresh rate is set to ',num2str(realrefrate),...
       'Hz while the proposed screen refresh rate is ',num2str(refrate),'Hz.']);
end

%Flicker1 frequency (Hz)
noisehz = str2num(answer{5});


%Flicker2 frequency (Hz)
flickerhz = str2num(answer{6});

%Number of Trials
trialnum = str2num(answer{7});

%% Code


%Colors
txtcolor = round([0 .6 .6]*255); %Teal
black = [0 0 0];
white = [255 255 255];
gray = 255*sqrt([.5 .5 .5]);
blackwhite{1} = black;
blackwhite{2} = white;

% Load fonts
myfont = '-bitstream-courier 14 pitch-bold-i-normal--0-0-0-0-m-0-ascii-0';
fontsize = 26;

%Define photocell placement
k = 0;
photorect = [0 0 100 90];
for m = 1:6
    pRect(m,:) = CenterRectOnPoint(photorect,50,50+k);
    k = k + 160;
end


%The following TRY, CATCH, END statement ends psychtoolbox if an error occurs
try
    %Open a fullscreen window on the first screen with black background
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');
    [wptr,windowRect] = PsychImaging('OpenWindow', whichScreen,gray);
    PsychGPUControl('FullScreenWindowDisablesCompositor', 1);
    
    %This vector defines the Flicker 1 frequency for our image
    noiseflic = [];
    for i=1:ceil(4*noisehz)
        noiseflic = [noiseflic 1 zeros(1,(round(refrate/noisehz)- 1))];
    end
    
    %This vector defines the Flicker 2 frequency for our image
    gaborflic = [];
    for i=1:ceil(4*flickerhz)
        gaborflic = [gaborflic 2*ones(1,round(refrate/2/flickerhz)) ones(1,round(refrate/2/flickerhz))];
    end
    
    %Inter-trial interval, 1500ms to 2000ms
    intertrial = 1.5 + rand(1,trialnum)*.5;
    
    %Calculate the number of frames in a cycle of an image flicker
    numCycleFrames = ceil(refrate*1.2) + ceil(refrate*rand(1,trialnum)*.8);
    
    Screen('TextFont',wptr,'Arial');
    Screen('TextSize',wptr,18);
    ShowCursor(0);	% arrow cursor
    sessiontext = 'Spacebar to continue';
    
    HideCursor;
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    Screen('DrawText',wptr, sessiontext,(xres - length(sessiontext)*9)/2,yres*(5/12),txtcolor);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    
    
    %Wait for spacebar
    FlushEvents('keyDown');
    [char,when] = GetChar; %Wait for keypress to continue
    notspace=1;
    while notspace
        switch char
            case ' '
                notspace =0;
            otherwise
                [char,when] = GetChar; %Wait for keypress to continue
                notspace =1;
        end
    end
    
    %Initialize timer
    tic;
    for trials = 1:trialnum
        
        %Display rush loops (Rush is apparently obsolete in PTB3, test this)
        Priority(MaxPriority(wptr)); %New to PTB3
        
        %Loop: Flicker 1 in photocells 1 and 5 from top
        %Flicker 2 in photocells 2 and 6 from top
        %Photocells 3 and 4 give one impulse for the beginning of a trial
        noisenum = 0;
        bwswitch = 1; 
        for i = 1:numCycleFrames(trials)
            if noiseflic(i)
                noisenum = noisenum + 1;
                bwswitch = mod(bwswitch,2) + 1; %Changes 1 to 2 and vica versa
            end
            if i == 1
                Screen(wptr,'FillRect',white,pRect([3 4],:)');
            else
                Screen(wptr,'FillRect',black,pRect([3 4],:)');
            end
            Screen(wptr,'FillRect',blackwhite{bwswitch},pRect([1 5],:)');
            Screen(wptr,'FillRect',blackwhite{gaborflic(i)},pRect([2 6],:)');
            Screen('Flip',wptr);
        end
        
        %Change all photocells to black after the trial has ended
        Screen(wptr,'FillRect',black,pRect');
        Screen('Flip',wptr);

        pause(intertrial(trials));
    end
        
catch me
    ShowCursor;
    Screen('CloseAll');
    rethrow(me); %rethrow reproduces the original error, stored in the object 'me'
end

ShowCursor;
Screen('CloseAll');
