function pdmexp3_public(varargin)
%% README

%We used this function to run a Perceptual Decision Making visual experiment for the
%two papers below. It gives an example of a steady-state visual evoked
%potential (SSVEP) experiment. See the two papers for a description of the
%experimental paradigm 

%THIS FUNCTION WILL NOT RUN WITHOUT SPECIFIC HARDWARE AND SOFTWARE
%
%We ran this on a Windows XP machine with Psychtoolbox 2. Subjects used a Cedrus
%button box.

%In order to replicate this experimental stimulus, you may need to port
%this code to Psychtoolbox 3 on Linux with a good video card, and then test 
%the frequency of the "photocell" squares in each corner of the screen to 
%ensure the entended frequency of the stimuli and thus possible SSVEP
%responses

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
%  12/?/2013     Michael Nunez                 Original code

%% Code
if nargin == 1
    rundemo = 1;
elseif nargin > 1
    error('Too many function inputs.');
else 
    rundemo = 0;
end
    
PTB2;

%Inputs Prompt and Output Setup
%Experimenter Prompt
Screenres = get(0,'Screensize');

prompt1={'Subject Number (must begin with letter):','Session Number:','Window Pointer:',...
    'Screen Length (x-axis):','Screen Width (y-axis):','Refresh Rate (fps):','Stimulus Frequencies (Hz <= Refresh Rate):',...
    'Noise Frequency:','Number of Trials:','Trials per Block:','Jitter Levels:','Noise Levels:','Cedrus Port [COM2 ...] or null'};
if rundemo
    def1={'SZZ_test','0','0',num2str(Screenres(3)),num2str(Screenres(4)),'120','15','8','36','90','[30 40 50]','[.3 .45 .6]',''};
    studytitle='PDM DEMO';
else
    def1={'SZZ_test','1','0',num2str(Screenres(3)),num2str(Screenres(4)),'120','15','8','540','90','[60 70 80]','[.3 .45 .6]',''};
    studytitle='PDM Experiment 3';
end

lineNo=1;
answer=inputdlg(prompt1,studytitle,lineNo,def1);
%Subject Number
subnum = answer{1};
%ExpSession Number
sesnum = str2num(answer{2});
output.sesnum = sesnum;
%Window Pointer / 'Home Screen'.  0 - the primary monitor; 1 - the secondary monitor.
whichScreen = str2num(answer{3});
%Screen resolution on the x-axis
xres = str2num(answer{4});
output.xres = xres;
%Screen resolution on the y-axis
yres = str2num(answer{5});
output.yres = yres;
%This should be the same as the Refresh Rate shown in the Display
%Properties on the computer.  Always check before running the experiment to
%match flicker frequency.
%This code is currently set up to only handle multiples of 60 fps.
refrate = str2num(answer{6});
realrefrate = Screen(0,'FrameRate');
if refrate ~= Screen(0,'FrameRate')
    error(['The real screen refresh rate is set to ',num2str(realrefrate),...
        'Hz while the proposed screen refresh rate is ',num2str(refrate),'Hz.']);
end
output.refrate = refrate;
%Stimulus frequencies(Hz)
flickfreqs = str2num(answer{7});
output.flickfreqs = flickfreqs;
if ~all(round(refrate./flickfreqs) == (refrate./flickfreqs))
    %error('The stimulus frequencies should be divisors of the refresh rate.');
end
%Noise frequency (Hz)
noisehz = str2num(answer{8});
output.noisehz = noisehz;
if round(refrate/noisehz) ~= refrate/noisehz
    error('The noise frequency should be divisor of the refresh rate.');
end
%Number of Trials
trialnum = str2num(answer{9});

%Trials per block
output.tperb = str2num(answer{10});
if output.tperb > trialnum
    output.tperb = trialnum;
end
block = 1;
%Jitter levels
jitterlvls = str2num(answer{11});
output.jitterlvls = jitterlvls;

%noise levels
noiselvls = str2num(answer{12});
output.noiselvls = noiselvls;

%Number of trials should be a multiple of the number of cells
ncells = length(jitterlvls)*length(noiselvls);
if round(trialnum/ncells) ~= trialnum/ncells
    error(['The number of trials should be a multiple of ',num2str(ncells)]);
end

%Cedrus Handle
cport = answer{13};
chandle = CedrusResponseBox('Open',cport);

%% Code
%Subject Prompt
prompt2={'What is your gender? (''f'' or ''m'')',...
    'Age:','Do you consider yourself right handed, left handed, or both?  (''r'',''l'', or''b'')',...
    'Do you have near 20/20 vision or is your vision corrected to near 20/20? (''y'' or ''n'')',...
    'Do you have any personal or family history of epilepsy? (''y'' or ''n'')'
    };
demographtitle='Subject Demographics';
lineNo=1;
subdemo=inputdlg(prompt2,demographtitle,lineNo);
switch subdemo{5}
    case 'n'
    otherwise
    error('You have indicated that you have a personal or family history of epilepsy. This experiment involves a fast flickering image. It is recommended that you NOT participate in this study due to a possible risk of seizure.  Please discuss your options with the experimenters.');
end
output.gender = subdemo{1};
output.age = str2num(subdemo{2});
output.hand = subdemo{3};
output.vision = subdemo{4};

%Get date and time that the session begins
output.date = date;
output.start_time = clock;
    
%number of rows and columns of image
nCols = 600;
nRows = 600;

%Initialize estimated accuracy vector, for speed
estcorrect = zeros(1,trialnum);

%Keyboard keypress variables
advancechar = ' ';
escapechar = 27;

%Lower bound to wait between trials, in seconds
lboundwait = 2.5;

% Initialize the sound driver:
InitializePsychSound;

%Load sounds
[ygood, Fsgood] = wavread('chimes.wav');
[ybad,  Fsbad ] = wavread('critical.wav');
ygood = ygood';
ybad = ybad'; 

% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
goodhand = PsychPortAudio('Open', [], [], 0, Fsgood, size(ygood,1));
badhand = PsychPortAudio('Open', [], [], 0, Fsbad, size(ybad,1));

% Fill the audio playback buffers with the audio data
PsychPortAudio('FillBuffer', goodhand , ygood);
PsychPortAudio('FillBuffer', badhand, ybad);

%Flush Cedrus Events
CedrusResponseBox('FlushEvents',chandle);


%The following TRY, CATCH, END statement ends psychtoolbox if an error
%occurs
try
    %opens home window
    [wptr windowRect]= Screen(whichScreen ,'OpenWindow',[255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)]);

    %sets size of gabor field that will be pasted onto Screen
    imageRect=SetRect(0,0,nCols,nRows);
    destRect=CenterRect(imageRect,windowRect);

    %Define photocell placement
    PhotoSize = 75;
    photorect = [0 0 PhotoSize PhotoSize];
    %Photo{4} = [0 0 PhotoSize PhotoSize];
    %Photo{3} = [xres-PhotoSize 0 xres PhotoSize];
    Photo{2} = [xres-PhotoSize yres-PhotoSize xres yres];
    Photo{1} = [0 yres-PhotoSize PhotoSize yres];
    
    %Creates windows for photocells
    owPCWhite = Screen(wptr, 'OpenOffScreenWindow', [255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)], [0 0 PhotoSize PhotoSize]);
    Screen(owPCWhite,'FillRect',[255 255 255],[0 0 PhotoSize PhotoSize]);
    owPCBlack = Screen(wptr, 'OpenOffScreenWindow', [255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)], [0 0 PhotoSize PhotoSize]);
    Screen(owPCBlack,'FillRect',[0 0 0],[0 0 PhotoSize PhotoSize]);

    %Create cross
    CrossWidth = 4;
    crossx1 = (round(nRows/2 - 3*CrossWidth/2):round(nRows/2 + 3*CrossWidth/2));
    crossy1 = (round((nCols)/2 - CrossWidth/2):round((nCols)/2 + CrossWidth/2));
    crossx2 = (round(nRows/2 - CrossWidth/2):round(nRows/2 + CrossWidth/2));
    crossy2 = (round((nCols)/2 - 3*CrossWidth/2):round((nCols)/2 + 3*CrossWidth/2));
    
    %Create white cross for flickering image
    whitecross = ones(600)*255*sqrt(.5); %Approximately accounts for monitor gamma
    whitecross(crossx1,crossy1) = ones(length(crossx1),length(crossy1))*255;
    whitecross(crossx2,crossy2) = ones(length(crossx2),length(crossy2))*255;
    
    %Creates window for white fixation cross with black photocells
    ow3 = Screen(wptr,'OpenOffScreenWindow',[255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)]);
    Screen(ow3,'PutImage',whitecross,destRect);
    Screen(ow3,'FillRect',[1 1 1],Photo{1});
    Screen(ow3,'FillRect',[1 1 1],Photo{2});
    %Screen(ow3,'FillRect',[1 1 1],Photo{3});
    %Screen(ow3,'FillRect',[1 1 1],Photo{4});

    %Black cross for beginning of stimulus
    blackcross = whitecross;
    blackcross(whitecross == 255) = 1;
    
    %Creates window for black fixation cross with black photocells
    ow2 = Screen(wptr,'OpenOffScreenWindow',[255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)]);
    Screen(ow2,'PutImage',blackcross,destRect);
    Screen(ow2,'FillRect',[1 1 1],Photo{1});
    Screen(ow2,'FillRect',[1 1 1],Photo{2});
    %Screen(ow2,'FillRect',[1 1 1],Photo{3});
    %Screen(ow2,'FillRect',[1 1 1],Photo{4});

    %Creates a window of a blank gray screen with black photocells
    ow4 = Screen(wptr,'OpenOffScreenWindow',[255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)]);
    Screen(ow4,'FillRect', 255*sqrt(.5), [0 0 xres yres]);
    Screen(ow4,'FillRect',[1 1 1],Photo{1});
    Screen(ow4,'FillRect',[1 1 1],Photo{2});
    %Screen(ow4,'FillRect',[1 1 1],Photo{3});
    %Screen(ow4,'FillRect',[1 1 1],Photo{4});
    
    %This matrix defines the flicker frequencies for our image
    stimflic = zeros( length(flickfreqs) , (round(.75*refrate)+ceil(refrate*2)) );
    for j=1:length(flickfreqs)
        particflick = [];
        for i=1:ceil(flickfreqs(j)*2)
            particflick = [particflick 2*ones(1,ceil(refrate/(2*flickfreqs(j)))) ones(1,floor(refrate/(2*flickfreqs(j))))];
        end
        stimflic(j,:) = [ones(1,round(.75*refrate)) particflick(1:ceil(refrate*2))]; %Does not display flickering bars until 750ms into stimulus 
    end
    
    %This vector defines the noise frequency for our image
    noiseflic = [];
    for i=1:ceil(noisehz*3)
        noiseflic = [noiseflic 1 zeros(1,(round(refrate/noisehz)- 1))];
    end
    noiseflic = noiseflic(1:ceil(refrate*3));
    
    %Set seed based on the time. Backwards compatible with older MATLAB
    %versions
    output.seed = sum(100*clock);
    rand('state',output.seed);
    
    %Set seed based on subject number.  Gives one the ability to replace
    %missing data
%     output.seed = str2num(subnum(2:end));
%     if isempty(output.seed)
%         output.seed = 0;
%         rng(0);
%     end
%     rng(output.seed);
    
    %Randomize frequency display
    whichfreq = [];
    for freq=1:length(flickfreqs)
        whichfreq = [whichfreq freq*ones(1,ceil(trialnum/length(flickfreqs)))];
    end
    whichfreq = whichfreq(randperm(length(whichfreq)));
    output.whichfreq = whichfreq;
    
    %Define jitter and noise vectors, ensure even cell counts
    jitnoise = [];
    for j=1:ncells
        jitnoise = [jitnoise ones(1,trialnum/ncells)*j];
    end
    jitnoise = jitnoise(randperm(length(jitnoise)));
    jittervec = zeros(1,length(jitnoise));
    noisevec = zeros(1,length(jitnoise));
    cellmat = reshape(1:ncells,length(jitterlvls),length(noiselvls));
    for j=1:length(jitterlvls)
        for k=1:length(noiselvls)
            jittervec(jitnoise == cellmat(j,k)) = jitterlvls(j);
            noisevec(jitnoise == cellmat(j,k)) = noiselvls(k);
        end
    end
    for j=1:length(noiselvls)
    end
    output.jittervec = jittervec;
    output.noisevec = noisevec;

    %Get vector of jitter values
    directions = [ones(1,round(trialnum/2)) ones(1,round(trialnum/2))*2];
    directions = directions(randperm(length(directions)));
    output.directions = directions;

    cut = 0; %Counter for ESC
    
    %Creates windows for fullimage
    for h=1:2
        for d=1:ceil(noisehz*3) %3 seconds of display
            owf{h,d} = Screen(wptr,'OpenOffScreenWindow',[255*sqrt(.5) 255*sqrt(.5) 255*sqrt(.5)]);
            %Place bar photocells in these offscreen windows
            if h==1
                Screen(owf{h,d},'FillRect',[1 1 1],Photo{2});
            else
                Screen(owf{h,d},'FillRect',[255 255 255],Photo{2});
            end
        end
    end
    
    %Set rush loop prioritys to Max
    priorityLevel=MaxPriority(whichScreen,'WaitBlanking');
    
    %Calculate the number of frames in a cycle of an image flicker
    numCycleFrames1 = round(.75*refrate) + (round(refrate)+round(refrate*rand(1,trialnum)));
    numCycleFrames2 = round(.5*refrate);

    %Output stimulus display time in milliseconds
    output.stimtime = (numCycleFrames1./60)*1000;
    
    %Rush Loop a: Black fixation cross for 750 ms
        rushloopa = {
            'Screen(''CopyWindow'', ow2, wptr);'
            'for frames = 1:round(.75*refrate);'
                'Screen(wptr,''WaitBlanking'');'
            'end;'
        };

        %Rush Loop 1: Noise for 750ms, noise with bars for 1000ms - 2000ms, accept responses
        rushloop1 = {
            'Screen(wptr,''WaitBlanking'');'
            'for i = 1:numCycleFrames1(trials);'
                'if noiseflic(i);'
                    'previousnoise = previousnoise + 1;'
                    'Screen(''CopyWindow'',owf{stimflic(whichfreq(trials),i),previousnoise},wptr);'
                    'Screen(''CopyWindow'',owPCWhite,wptr,photorect,Photo{1});'
                'else;'
                    'Screen(''CopyWindow'',owf{stimflic(whichfreq(trials),i),previousnoise},wptr);'
                    'Screen(''CopyWindow'',owPCBlack,wptr,photorect,Photo{1});'
                'end;'
                'Screen(wptr,''WaitBlanking'');'
            'end;'
        };

        %Rush Loop 2: Keep displaying black fixation cross (only) for 250ms
        rushloop2 = {
            'Screen(''CopyWindow'', ow2, wptr);'
            'for frames = 1:round(refrate/4);'
                'Screen(wptr,''WaitBlanking'');'
            'end;'
        };
    
    Screen(wptr,'TextFont','Arial');
    Screen(wptr,'TextSize',18);
    ShowCursor(0);	% arrow cursor
    sessiontext = 'Loading images...';
    sessiontext2 = 'The experiment will begin shortly';
    sessiontext3 = 'The experiment has started! Good luck!';
    
    Screen('CopyWindow', ow4, wptr);
    HideCursor;
    Screen(wptr, 'DrawText',sessiontext,(xres - length(sessiontext)*9)/2,yres/2,[1 1 1]);
    
    
    %Load first block's bar field images
    for b=1:output.tperb
        barimage(:,:,b) = makebardirec(600,6,directions(b),jittervec(b),0); %Get the gabor image color data using makebarimage.m
        barimage(barimage == .5) = .15;
    end
    
    %Display second text screen
    Screen('CopyWindow', ow4, wptr);
    Screen(wptr, 'DrawText',sessiontext2,(xres - length(sessiontext2)*9)/2,yres/2,[1 1 1]);
    
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
    
    %Display third text screen
    Screen('CopyWindow', ow4, wptr);
    Screen(wptr, 'DrawText',sessiontext3,(xres - length(sessiontext3)*9)/2,yres/2,[1 1 1]);
    
    %Initialize timer
    tic;

    for trials = 1:trialnum
      if ~cut %ESC key track
        
        imageset = cell(2,ceil(noisehz));
        trialind = trials-(block-1)*output.tperb;
        for f=1:ceil(noisehz)
            %Create image of 5% contrast (defined as contrast ratio
            %whitest/darkest)
            imageset{1,f} = (.5 - (noisevec(trials)/2)) + noisevec(trials)*rand(600);
            imageset{2,f} = barimage(:,:,trialind) + imageset{1,f};
            %transformation for barimage color values, takes into account
            %monitor gamma
            imageset{1,f} = 255*sqrt(imageset{1,f}); %The square root is in order to account for monitor gamma. That is, the monitor approximately squares the input stimulus color value
            imageset{2,f} = 255*sqrt(imageset{2,f});
            %Create black cross ontop of fullimage
            imageset{1,f}(crossx1,crossy1) = ones(length(crossx1),length(crossy1));
            imageset{1,f}(crossx2,crossy2) = ones(length(crossx2),length(crossy2));
            imageset{2,f}(crossx1,crossy1) = ones(length(crossx1),length(crossy1));
            imageset{2,f}(crossx2,crossy2) = ones(length(crossx2),length(crossy2));
        end
        
        %Creates windows for fullimage
        for q=1:3 %3 seconds of frames
            for d=1:ceil(noisehz)
                frameind = d + noisehz*(q-1);
                Screen(owf{1,frameind},'PutImage',imageset{1,d},destRect);
                Screen(owf{2,frameind},'PutImage',imageset{2,d},destRect);
            end
        end
        
        %Initialize previous noise
        previousnoise = 0;
        
        %Wait at least lboundwait seconds between trials
        output.elapsedtime(trials) = toc;
        if output.elapsedtime(trials) < lboundwait
            pause(lboundwait-output.elapsedtime(trials));
        end
        output.fixedtime(trials) = toc;
        
        CedrusResponseBox('FlushEvents',chandle);
        
        %Display rush loops
        Rush(rushloopa,priorityLevel);
        Rush(rushloop1,priorityLevel);
        Rush(rushloop2,priorityLevel);

        %Timer to calculate time between the last trial and the next
        tic;

        %Show the blank gray Screen
        Screen('CopyWindow', ow4, wptr);
        Screen(wptr,'WaitBlanking');
        
        %Play feedback sound
        evt = CedrusResponseBox('GetButtons',chandle);
        if isempty(evt)
            correct = 0;
        elseif (evt.button == 5 && directions(trials) == 1) || ...
                (evt.button == 3 && directions(trials) == 2)
            correct = 1;
        else
            correct = 0;
        end
        estcorrect(trials) = correct;
        if correct
            %sound(ygood, Fsgood);
            PsychPortAudio('Start', goodhand , 1, 0, 1);
        else
            %sound(ybad, Fsbad);
            PsychPortAudio('Start', badhand , 1, 0, 1);
        end

        if ~cut
            if trials == trialnum
                %Show ending screen for 5 seconds
                percorrect = sum(estcorrect((trials-output.tperb+1):trials))/output.tperb;
                endtext = ['Done!  ',...
                    num2str(round(percorrect*100)),'% correct this block. Thank you for participating!'];
                Screen('CopyWindow', ow4, wptr,[0 0 xres yres]);
                Screen(wptr, 'DrawText',endtext,(xres - length(endtext)*9)/2,yres/2,[1 1 1]);
                %Wait for spacebar to end program
                FlushEvents('keyDown');
                [char,~] = GetChar; %Wait for keypress to continue
                notspace=1;
                while notspace
                    switch char
                        case advancechar
                            notspace =0;
                        otherwise
                            [char,~] = GetChar; %Wait for keypress to continue
                            notspace =1;
                    end
                end
            elseif trials/output.tperb == round(trials/output.tperb)
                 %Take a break every 'output.tperb' trials and show ending Screens
                percorrect = sum(estcorrect((trials-output.tperb+1):trials))/output.tperb;
                trialtext = ['Block ',num2str(block),' complete!  ',...
                    num2str(round(percorrect*100)),'% correct this block. You may now take a break!'];
                block = block + 1;
                trialtext2 = 'Please wait for the experimenter';
                Screen('CopyWindow', ow4, wptr,[0 0 xres yres]);
                Screen(wptr, 'DrawText',trialtext,(xres - length(trialtext)*9)/2,yres/2,[1 1 1]);
                Screen(wptr, 'DrawText',trialtext2,(xres - length(trialtext2)*9)/2,yres/2 + 32,[1 1 1]);
                
                %Load first block's images
                for b=1:output.tperb %Load first block's images
                    barimage(:,:,b) = makebardirec(600,6,directions(b+trials),jittervec(b+trials),1); %Get the gabor image color data using makebarimage.m
                    barimage(barimage == .5) = .15;
                end
                
                %Wait for spacebar
                FlushEvents('keyDown');
                [char,~] = GetChar; %Wait for keypress to continue
                notspace=1;
                while notspace
                    switch char
                        case advancechar
                            notspace =0;
                            Screen('CopyWindow', ow4, wptr);
                            Screen(wptr, 'DrawText',sessiontext3,(xres - length(sessiontext3)*9)/2,yres/2,[1 1 1]);
                            %Timer to calculate time between the last trial and the next
                            tic;
                        case escapechar %Escape from experiment and save current data (for experimenter)
                            notspace =0;
                            RestoreScreen(whichScreen);
                            ShowCursor;
                            Screen('Closeall');
                            output.ESC_time = clock;
                            output.estcorrect = estcorrect;
                            eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
                            if ~exist('data','dir')
                                mkdir('data');
                            end
                            eval(['save(''data/',subnum,'_Exp_',num2str(output.sesnum),'_',date,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
                            warning on all;
                            CedrusResponseBox('Close',chandle);
                            PsychPortAudio('Close', goodhand);
                            PsychPortAudio('Close', badhand);
                            return
                        otherwise
                            [char,when] = GetChar; %Wait for keypress to continue
                            notspace =1;
                    end
                end
            end
            
        end
      end  
    end
catch me
    RestoreScreen(whichScreen);
    ShowCursor;
    Screen('Closeall');
    output.error_time = clock;
    output.estcorrect = estcorrect;
    eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
    if ~exist('data','dir')
        mkdir('data');
    end
    eval(['save(''data/',subnum,'_Exp_',num2str(output.sesnum),'_',date,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
    CedrusResponseBox('Close',chandle);
    PsychPortAudio('Close', goodhand);
    PsychPortAudio('Close', badhand);
    rethrow(me); %rethrow reproduces the original error, stored in the object 'me'
end

RestoreScreen(whichScreen);
ShowCursor;
Screen('Closeall');

%Output time finished
output.finish_time = clock;

%Estimated accuracy
output.estcorrect = estcorrect;

eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
if ~exist('data','dir')
    mkdir('data');
end
    eval(['save(''data/',subnum,'_Exp_',num2str(output.sesnum),'_',date,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
warning on all;
CedrusResponseBox('Close',chandle);
%PsychPortAudio('Close', goodhand);
%PsychPortAudio('Close', badhand);


