function mondichNRexpt

% normalization reweighting experiment with mon and dichoptic presentation
% code based on surround suppression experiment from 2014
% this version for ViewPixx, 36 pix/cm dot pitch
% four conditions - 5Hz component, 7Hz component, both mon, dichoptic
% DHB 1/3/20

clear
close all;

pID = inputdlg('Enter participant ID number');
E.subj = pID{1,1};
useVP = 0;                  % hardware stereo doesn't seem to work in M16 mode
atYNIC = 0;
blockstorun = 3;

E.nconds = 4;                           % 5Hz component, 7Hz component, mon combination, dich combination
E.ntrials = 16;                         % total number of trials per condition in each block
E.ntrialsperblock = E.nconds*E.ntrials;
E.totalblocks = 3;
E.totalreps = E.ntrials*E.totalblocks;

E.levelsC5 = [48 0 48 48]./100;
E.levelsC7 = [0 48 48 48]./100;

W = what; E.exptpath = strcat(W.path,'/'); clear W;
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 1);

if atYNIC
    ST.npixelsperdegree = 44;       % MEG projector pixel resolution at 85cm
    addpath(strcat(E.exptpath,'/ppdev-mex-master'));
else
    ST.npixelsperdegree = 36;       % VPixx pixel resolution at 57cm
end
ST.gratingsize = ST.npixelsperdegree*2;         %
ST.gratingradiusfactor = ST.npixelsperdegree*2;  % previously 128
ST.SF = 2;
ST.ncycles = ST.SF*ST.gratingsize/ST.npixelsperdegree;

ST.duration = 6;
ST.ITI = 3;         % always leave n seconds between trials
ST.TFtarget = 5;
ST.TFmask = 7;

if (~exist(strcat(E.exptpath, 'Results/'),'dir'))
    mkdir(strcat(E.exptpath, 'Results/'));
end
fname = strcat(E.exptpath, 'Results/', E.subj, 'trialorderNR.mat');
if exist(fname)
    load(fname);
    save(strcat(E.exptpath, 'Results/', E.subj, 'trialorderNRold.mat'),'R','E');
else
    R.subj = E.subj;
    R.resps = zeros(E.nconds,E.totalreps);
    R.tno = 0.*(1:E.nconds);
    
    for n = 1:E.totalblocks
        R.trialorder(n,:) = randperm(E.ntrialsperblock);
        R.targetleft(n,:) = rem(1:E.ntrialsperblock,2);
        R.orientations(n,:,:) = 360.*(rand(36,E.ntrialsperblock));
    end
    R.currentrep = 0;
    
    save(fname, 'R', 'E');
end

WaitSecs(0.01);                % important to load in the MEX file before the expt starts
GetSecs;
InitializePsychSound;
tr = PsychPortAudio('Open',[],[],[],[],1);
PsychPortAudio('FillBuffer', tr, MakeBeep(440*sqrt(2),0.05,44100).*0.5);

    
try                 % start the 'try/catch' loop
    
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    PsychGPUControl('SetDitheringEnabled', 0);
    
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    if useVP        % using a ViewPixx or ProPixx
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        if ~atYNIC
            Datapixx('DisableVideoScanningBacklight');      % optionally, turn it off first, in case the refresh rate has changed since startup
            Datapixx('EnableVideoScanningBacklight');       % Only required if a VIEWPixx.
            Datapixx('EnableVideoStereoBlueline');
            Datapixx('SetVideoStereoVesaWaveform', 2);      % If driving NVIDIA glasses
        end
        PsychImaging('AddTask', 'General', 'EnableDataPixxM16Output');
        if Datapixx('IsViewpixx3D')
            Datapixx('EnableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        
        if atYNIC   % open parallel port
            ppdev_mex('Open', 1);
            ppdev_mex('Write', 1, 0);
        end

        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
        [w, winRect] = PsychImaging('OpenWindow', screenNumber, 0, [], [], [], 0, 0, kPsychNeedFastBackingStore);
        Screen('LoadNormalizedGammaTable', w, linspace(0,1,256)'*ones(1,3), 0); % THIS IS THE IMPORTANT THING TO DO, NOTE THE LAST ARGUMENT IS 0.

        HideCursor;
        ST.greylevel = 0.5;
        if ~atYNIC
            ST.greylevel = doimagegamma(0.5);
        end
    else
        rect = [1 1 1024 1024];
        [w, winRect] = Screen('OpenWindow',screenNumber,128,rect);
        ST.greylevel = round(doimagegamma(0.5)*255);
    end
    
    
    ST.blueRectLeftOn   = [0, winRect(4)-1, winRect(3)/4, winRect(4)];
    ST.blueRectLeftOff  = [winRect(3)/4, winRect(4)-1, winRect(3), winRect(4)];
    ST.blueRectRightOn  = [0, winRect(4)-1, winRect(3)*3/4, winRect(4)];
    ST.blueRectRightOff = [winRect(3)*3/4, winRect(4)-1, winRect(3), winRect(4)];
    ST.greyRect = winRect + [0 1 0 -1];
    ST.upperblackrect = [0, winRect(2), winRect(3), winRect(2)+1];
    ST.lowerblackrect = [0, winRect(4)-1, winRect(3), winRect(4)];
    
    Screen('FillRect', w, ST.greylevel);
    Screen('Flip', w);
    
    [width, height] = Screen('WindowSize', w);
    ifi=Screen('GetFlipInterval', w);
    ifims = ifi*1000*2;        % doubled because of frame interleaving
    
    ndots = 250;
    doffset = 320;
    dotangle = rand(1,ndots).*2.*pi;
    dotradius = rand(1,ndots).*ST.npixelsperdegree/4;
    [ST.dotcoords(1,:) ST.dotcoords(2,:)] = pol2cart(dotangle, dotradius);
    ST.dotcoords(1,51:100) = ST.dotcoords(1,51:100)+doffset;
    ST.dotcoords(1,101:150) = ST.dotcoords(1,101:150)-doffset;
    ST.dotcoords(1,151:200) = ST.dotcoords(1,151:200)-doffset;
    ST.dotcoords(1,201:250) = ST.dotcoords(1,201:250)+doffset;
    ST.dotcoords(2,51:100) = ST.dotcoords(2,51:100)+doffset;
    ST.dotcoords(2,101:150) = ST.dotcoords(2,101:150)+doffset;
    ST.dotcoords(2,151:200) = ST.dotcoords(2,151:200)-doffset;
    ST.dotcoords(2,201:250) = ST.dotcoords(2,201:250)-doffset;
    
    ST.dotlevels = rand(3,ndots);
    
    ST.dotsize = 8;
    ST.centre = [width/2, height/2];
    
    Screen('FillRect',w, 0.1, ST.greyRect);
    Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots
    Screen('Flip', w);
    
    ST.nframes = round(1000/ifims);       % this is 60Hz
    
    nframes = 60;          % one full second
    targetwaveform = -cos(2 .* ST.TFtarget .* (ifims:ifims:(nframes*ifims)) .* pi./1000);
    targetwaveform = (targetwaveform + 1)./2;
    maskwaveform = -cos(2 .* ST.TFmask .* (ifims:ifims:(nframes*ifims)) .* pi./1000);
    maskwaveform = (maskwaveform + 1)./2;
    
    trigwaveform1 = -cos(2 .* ST.TFtarget .* ((1/nframes):(1/nframes):(2*ST.duration.*nframes*(1/nframes))) .* pi);
    trigwaveform1 = round(48.*(trigwaveform1 + 1)./2);
    trigwaveform2 = -cos(2 .* ST.TFmask .* ((1/nframes):(1/nframes):(2*ST.duration.*nframes*(1/nframes))) .* pi);
    trigwaveform2 = round(48.*(trigwaveform2 + 1)./2);
    trigwaveform3 = trigwaveform1 + trigwaveform2;
    
    gratingangle = [(0:90:330) (45:90:330) (0:90:330) (22.5:45:360) (0:45:315) (22.5:45:360)];
    gratingangle = pi.*gratingangle./180;
    gratingradius(1:4) = ST.gratingradiusfactor;
    gratingradius(5:8) = ST.gratingradiusfactor.*2;
    gratingradius(9:12) = ST.gratingradiusfactor.*2.2;
    gratingradius(13:20) = ST.gratingradiusfactor.*3;
    gratingradius(21:28) = ST.gratingradiusfactor.*3.6;
    gratingradius(29:36) = ST.gratingradiusfactor.*4.4;
    
    [gratingcoords(1,:) gratingcoords(2,:)] = pol2cart(gratingangle, gratingradius);
    
    r1 = [1 1 ST.gratingsize ST.gratingsize];
    destRect = CenterRectOnPoint(r1, width*0.5, height*0.5);
    
    for n = 1:length(gratingcoords)
        rectlist(:,n) = CenterRectOnPoint(r1, width*0.5+gratingcoords(1,n), height*0.5+gratingcoords(2,n));
    end
    
    sgrating = mkgrating(ST.gratingsize, ST.ncycles, 0, 90, 1);
    window = make_soft_window(ST.gratingsize,ST.gratingsize,0.9);
    
    for n = 1:nframes              % target only
        comp = targetwaveform(n).*sgrating.*window.*E.levelsC5(1);
        comp = (1+comp)/2;
        if ~atYNIC
            comp = doimagegamma(comp);
        end
        targettexture(n) = Screen('MakeTexture', w, comp, [], [], 2);
        
        comp = maskwaveform(n).*sgrating.*window.*E.levelsC7(2);
        comp = (1+comp)/2;
        if ~atYNIC
            comp = doimagegamma(comp);
        end
        masktexture(n) = Screen('MakeTexture', w, comp, [], [], 2);
        
        comp = targetwaveform(n).*sgrating.*window.*E.levelsC5(1) + maskwaveform(n).*rot90(sgrating).*window.*E.levelsC7(2);
        comp = (1+comp)/2;
        if ~atYNIC
            comp = doimagegamma(comp);
        end
        plaidtexture(n) = Screen('MakeTexture', w, comp, [], [], 2);
    end
    
    comp = sgrating.*0;
    comp = (1+comp)/2;
    if ~atYNIC
        comp = doimagegamma(comp);
    end
    nulltexture = Screen('MakeTexture', w, comp, [], [], 2);
    
        PsychPortAudio('Start', tr);
    WaitSecs(0.5);
    
    % shows null stimuli to avoid flicker on initial frame of first trial
                Screen('FillRect',w, ST.greylevel, ST.greyRect);
                Screen('DrawTextures', w, nulltexture, [], rectlist);
                if useVP
                    Screen('FillRect', w, [255 255 255], ST.blueRectLeftOn);
                    Screen('FillRect', w, [0 0 0], ST.blueRectLeftOff);
                end
                Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w);
                lastflip = Screen('Flip', w);

                Screen('FillRect',w, ST.greylevel, ST.greyRect);
                Screen('DrawTextures', w, nulltexture, [], rectlist);
                if useVP
                    Screen('FillRect', w, [255 255 255], ST.blueRectRightOn);
                    Screen('FillRect', w, [0 0 0], ST.blueRectRightOff);
                end
                Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots
                Screen('DrawingFinished', w);
                lastflip = Screen('Flip', w, lastflip+ifi*0.5);
 
                
    alltriggertimes = [];
    
    breakcode = 0;
    currentblock = 0;
    
    while currentblock < blockstorun
        
        currentblock = currentblock + 1;
        R.currentrep = R.currentrep + 1;
        
        R.reportedifi(R.currentrep) = ifi;
                
        PsychPortAudio('Start', tr);
        
        exitcode = 0;
        while exitcode==0
            
            if useVP
                lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
            else
                lastflip = Screen('Flip', w);
            end
            
            [x,y,buttons] = GetMouse;
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if sum(buttons)>0 || keyIsDown
                exitcode = 1;
            end
        end
        
        
        for n = 1:3                         % 3 trigger pulses to indicate start of block
            if useVP
                if atYNIC
                    ppdev_mex('Write', 1, 1);
                else
                    Datapixx('SetDoutValues', transformindex(1));
                    Datapixx('RegWrRd');
                end
            end
            for nwaitframepairs = 1:3
                if useVP
                    lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
                else
                    lastflip = Screen('Flip', w);
                end
            end
            if useVP
                if atYNIC
                    ppdev_mex('Write', 1, 0);
                else
                    Datapixx('SetDoutValues', 0);
                    Datapixx('RegWrRd');
                end
            end
            for nwaitframepairs = 1:3
                if useVP
                    lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
                else
                    lastflip = Screen('Flip', w);
                end
            end
            
        end
        
        if useVP
            lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
        else
            lastflip = Screen('Flip', w);
        end
        trialoffset = lastflip;     % last trial was an infinitely long time ago
        trial = 0;
        
        while trial < E.ntrialsperblock

            trial = trial + 1;
            thistrial = R.trialorder(R.currentrep,trial);
            E.cond = ceil(thistrial/E.ntrials);
            isleft = R.targetleft(R.currentrep,trial);
            stimorient = squeeze(R.orientations(R.currentrep,:,trial));
            R.tno(E.cond) = R.tno(E.cond) + 1;

            switch E.cond
                case 1      % mon target only
                    stimlist1 = targettexture;
                    stimlist2 = nulltexture;
                    triglist1 = trigwaveform1;
                    triglist2 = trigwaveform1.*0;
                case 2      % mon  mask
                    stimlist1 = masktexture;
                    stimlist2 = nulltexture;
                    triglist1 = trigwaveform2;
                    triglist2 = trigwaveform2.*0;
                case 3      % mon target plus orthogonal mask
                    stimlist1 = plaidtexture;
                    stimlist2 = nulltexture;
                    triglist1 = trigwaveform3;
                    triglist2 = trigwaveform3.*0;
                case 4      % mon target plus dich mask
                    stimlist1 = targettexture;
                    stimlist2 = masktexture;
                    triglist1 = trigwaveform1;
                    triglist2 = trigwaveform2;
            end

            if ~isleft
                stimlistL = stimlist1;
                stimlistR = stimlist2;
                triglistL = triglist1;
                triglistR = triglist2;
            else
                stimlistL = stimlist2;
                stimlistR = stimlist1;
                triglistL = triglist2;
                triglistR = triglist1;
            end

            nframesL = length(stimlistL);
            nframesR = length(stimlistR);
            frametoswitch = round(2*rand*(ST.nframes*ST.duration));     %switch the fixation marker with a 50% probability on each trial

            while lastflip<(trialoffset+ST.ITI)         % leave a gap in between each trial
                if useVP
                    lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
                else
                    Screen('FillRect', w, ST.greylevel);
                    lastflip = Screen('Flip', w);
                end
            end
            if useVP
                if atYNIC
                    ppdev_mex('Write', 1, 200+E.cond*10+isleft);
                else
                    Datapixx('SetDoutValues', transformindex(200+E.cond*10+isleft));     % this trigger contains condition code
                    Datapixx('RegWrRd');
                end
                lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
                if atYNIC
                    ppdev_mex('Write', 1, 0);
                else
                    Datapixx('SetDoutValues', 0);
                    Datapixx('RegWrRd');
                end
                lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
            end
            
            starttime = lastflip;
            
            n = 0;
            while lastflip < (starttime + ST.duration)
                % this code to compensate for skipped frames
                n = round((lastflip - starttime)/(2*ifi)) + 1;
                frameindexL = mod(n-1,nframesL)+1;
                frameindexR = mod(n-1,nframesR)+1;
                frameindex = mod(n-1,ST.nframes)+1;
  
                Screen('FillRect',w, ST.greylevel, ST.greyRect);
                Screen('DrawTextures', w, stimlistL(frameindexL), [], rectlist, stimorient-45);
                if useVP
                    Screen('FillRect', w, [255 255 255], ST.blueRectLeftOn);
                    Screen('FillRect', w, [0 0 0], ST.blueRectLeftOff);
                end

                Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots

                if useVP
                    if atYNIC
                        ppdev_mex('Write', 1, triglistL(n));
                    else
                        Datapixx('SetDoutValues', transformindex(triglistL(n)));     % first trigger also contains condition code
                        Datapixx('RegWrRd');
                    end
                end

                Screen('DrawingFinished', w);
                lastflip = Screen('Flip', w, lastflip+ifi*0.5);

                Screen('FillRect',w, ST.greylevel, ST.greyRect);
                Screen('DrawTextures', w, stimlistR(frameindexR), [], rectlist, stimorient+45);
                if useVP
                    Screen('FillRect', w, [255 255 255], ST.blueRectRightOn);
                    Screen('FillRect', w, [0 0 0], ST.blueRectRightOff);
                end
                Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots

                if useVP
                    if atYNIC
                        ppdev_mex('Write', 1, triglistR(n));
                    else
                        Datapixx('SetDoutValues', transformindex(triglistR(n)));     % first trigger also contains condition code
                        Datapixx('RegWrRd');
                    end
                end
                Screen('DrawingFinished', w);
              
                lastflip = Screen('Flip', w, lastflip+ifi*0.5);
                frametime(n) = lastflip - starttime;
                               
                if n==frametoswitch
                    dotangle = rand(1,ndots).*2.*pi;
                    dotradius = rand(1,ndots).*ST.npixelsperdegree/4;
                    [ST.dotcoords(1,:) ST.dotcoords(2,:)] = pol2cart(dotangle, dotradius);
                    ST.dotcoords(1,51:100) = ST.dotcoords(1,51:100)+doffset;
                    ST.dotcoords(1,101:150) = ST.dotcoords(1,101:150)-doffset;
                    ST.dotcoords(1,151:200) = ST.dotcoords(1,151:200)-doffset;
                    ST.dotcoords(1,201:250) = ST.dotcoords(1,201:250)+doffset;
                    ST.dotcoords(2,51:100) = ST.dotcoords(2,51:100)+doffset;
                    ST.dotcoords(2,101:150) = ST.dotcoords(2,101:150)+doffset;
                    ST.dotcoords(2,151:200) = ST.dotcoords(2,151:200)-doffset;
                    ST.dotcoords(2,201:250) = ST.dotcoords(2,201:250)-doffset;
                end
            
            end
            
            allframetimes(trial,1:length(frametime)) = frametime;
            trialoffset = lastflip;
            
            if useVP
                if atYNIC
                    ppdev_mex('Write', 1, 0);
                else
                    Datapixx('SetDoutValues', 0);
                    Datapixx('RegWrRd');
                end
                for nwaitframepairs = 1:3
                    lastflip = dogoggleflippair(w,ST,ifi,w,lastflip);
                end
            else
                Screen('FillRect', w, ST.greylevel);
                lastflip = Screen('Flip', w);
            end
            
            [keyIsDown, secs, keyCode] = KbCheck;
            
            if keyCode(KbName('RightArrow'))
                breakcode = 1;
                trial = 1000;
                E.block = 100;
                currentblock = 100;
            end
            
        end
        
        if ~breakcode
            save(fname, 'R', 'E');       % save every block of trials
            save(strcat(E.exptpath, 'Results/', E.subj,'TimesBlock',num2str(R.currentrep),'.mat'), 'allframetimes', 'alltriggertimes');
        end
        
        endofblock = lastflip;
        while (lastflip < (endofblock+2))
            
            if useVP
                
                Screen('FillRect', w, ST.greylevel, ST.greyRect);
                Screen('FillRect', w, [255 255 255], ST.blueRectLeftOn);
                Screen('FillRect', w, [0 0 0], ST.blueRectLeftOff);
                
                Screen('DrawText', w, 'End of block', width/2, height/2);
                
                Screen('DrawingFinished', w);
                lastflip = Screen('Flip', w, lastflip+ifi*0.5);
                
                Screen('FillRect',w, ST.greylevel, ST.greyRect);
                Screen('FillRect', w, [255 255 255], ST.blueRectRightOn);
                Screen('FillRect', w, [0 0 0], ST.blueRectRightOff);
                
                Screen('DrawText', w, 'End of block', width/2, height/2);
                
                Screen('DrawingFinished', w);
                lastflip = Screen('Flip', w, lastflip+ifi*0.5);
                
            else
                Screen('FillRect', w, ST.greylevel);
                Screen('DrawText', w, 'End of block', width/2, height/2);
                lastflip = Screen('Flip', w);
            end
            
        end
        
        
    end
    
catch
    
    lasterr
    
end

Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

Screen('Flip', w);
Screen('Close',nulltexture);
Screen('Close',targettexture(:));
Screen('Close',masktexture(:));
Screen('Close',plaidtexture(:));
ShowCursor;

if useVP        % using a ViewPixx
    if ~atYNIC
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
    end
end
Screen('CloseAll');
if useVP
    Datapixx('Close');
end

PsychPortAudio('Close', tr);

if atYNIC   % close parallel port
    ppdev_mex('Close', 1);
end

s = size(allframetimes);
for a = 1:s(1)
    temp = allframetimes(a,:);
    framelengths(a,:) = temp(2:end) - temp(1:end-1);
end

a = find(framelengths(:)>(1.1*2*ifi));
percentskippedframes = 100*length(a)./length(framelengths(:))

end
%--------------------------------------------------------------------------
function mouseloop

exitcode = 0;

while exitcode==0
    [x,y,buttons] = GetMouse;
    
    if sum(buttons)>0
        exitcode = 1;
    end
end

end
%--------------------------------------------------------------------------------------------------
function imag1 = mkgrating(Regionsize, f, o, p, c)

%TSM; 26.6.03
% modified by DHB to make single component gratings only, scaled from -1 to 1
% f is spatial frequency, scaled as cycles per image
% o is orientation (degrees), p is phase (degrees relative to centre), c is contrast

p = p*pi/180;
o = o*2*pi/360;		% convert from degrees to radians
f = f/Regionsize;
x0 = ((Regionsize+1)/2);
y0 = x0;

u = f .* cos(o) * 2 * pi;
v = f .* sin(o) * 2 * pi;

imag1 = zeros(Regionsize, Regionsize);
[xx, yy] = meshgrid(1:Regionsize, 1:Regionsize);

imag1(:,:) = (c .* sin(u .*(xx-x0) + v.*(yy-y0) + p));

end
%--------------------------------------------------------------------------
function dB = d_PercentTodB(perc)

%converts from % contrast to contrast in dB

dB = 20 * log10(perc);

end
%--------------------------------------------------------------------------------------------------
function [perc] = d_dBtoPercent(dB)

%converts from dB to % contrast

perc = 10.^(dB/20);

end
%--------------------------------------------------------------------------
function kbloop

exitcode = 0;

while exitcode==0
    [keyIsDown, secs, keyCode] = KbCheck;       % also monitor keyboard for breaks
    
    if sum(keyCode(79:80))>0        % respond to left and right arrows
        exitcode = 1;
    end
end

end
%--------------------------------------------------------------------------------------------------
function mask = make_soft_window(W,H,D)

% Mark's code for making a raised cosine window

% SYNTAX: mask = make_soft_window(W,H,[D])
% ends an array 'mask' that is 1 inside the circular window, shading to zero outside
% W, H are the width and height of the whole (rectangular or square) array, in pixels
% Diameter of the soft window at half-height defaults to 0.90 units
%    where 1 unit = image width or height (whichever is the smaller)
% Smoothing is by convolution with a cosine half-cycle of width 0.1 units
% Optional parameter D specifies this diameter (in relative units, range 0 -> 1)
% MAG, 27.2.04

%soft window parameters
if nargin<3, D = 0.9; end % sets default diameter to 0.9
radius = min(W*D/2,H*D/2);% radius in pixels
blur = 2*(min(W/2,H/2) - radius);  % blur half-cycle
L = blur;
X1 = [-L/2:L/2];

% 1-D blur function (applied twice, in x and y)
WinKernel = cos(X1*pi/L); % half-cycle cosine
%image coordinates - X and Y arrays
X = [1:W] - W/2;
Y = [1:H] - H/2;
xx = repmat(X,H,1);
yy = repmat(Y',1,W);

% make circular soft window
mask = single((xx.*xx + yy.*yy) < radius^2); % logical 0 outside the circle,1 inside it
mask = conv2(WinKernel,WinKernel,mask,'same'); 	% smooth the mask
mask = mask/max(max(mask));						% scale the mask 0-1
% figure(2);plot(X,mask(H/2,:),'r-',Y,mask(:,W/2));
mask = double(mask);
end
%--------------------------------------------------------------------------------------------------
function output = doimagegamma(i)

% gamma corrects the stimuli before sending to Bits++
% adapted from Mark's code, DHB 29.01.08

% parameters from last gamma correct, 12/8/14 on VPixx in M16 mode
k = 0.5344;
Lmax = 101.8721;
j0 = 0.1292;
gamma = 1.9252;
%%%%%

i0 = 0;
imax = 1;                                               % Bits++ always scaled between 0 and 1
imean = (i0+imax)/2;
jmax = 1;

Lmin = k + (Lmax-k)*(max(-j0,0)/(jmax-j0) ).^gamma;     % Eqn 2, with j set to 0, to get Lmin
Lmin = max(Lmin,0);                                     % ensure Lmin not <0
Lmean = (Lmin+Lmax)/2;
L = Lmean + (Lmax-Lmean)*(i-imean)/(imax-imean);        % desired luminance values Eqn 4
j = ((L - k)/(Lmax-k)).^(1/gamma)*(jmax - j0) + j0;     % These are the gamma-corrected lut values, j: Eqn 3
output = max(j,j0);                                     % Eqn 3 conditional
output = double(output);

end
%--------------------------------------------------------------------------
function lastflip = dogoggleflippair(w,ST,ifi,overlay,lastflip)

Screen('FillRect', w, ST.greylevel, ST.greyRect);
Screen('FillRect', overlay, [255 255 255], ST.blueRectLeftOn);
Screen('FillRect', overlay, [0 0 0], ST.blueRectLeftOff);

Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots

Screen('DrawingFinished', w);
lastflip = Screen('Flip', w, lastflip+ifi*0.5);

Screen('FillRect',w, ST.greylevel, ST.greyRect);
Screen('FillRect', overlay, [255 255 255], ST.blueRectRightOn);
Screen('FillRect', overlay, [0 0 0], ST.blueRectRightOff);
Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots

Screen('DrawingFinished', w);
lastflip = Screen('Flip', w, lastflip+ifi*0.5);

end
%--------------------------------------------------------------------------
function output = transformindex(input)

% fixes the binary inputs for the EEG amplifier because the pins are in a different order from the ViewPixx
% desired numbers must be <256
% DHB 18/8/14

truebits = 2.^(2:2:24);
dn = dec2bin(input,length(truebits));
output = 0;
for m = 1:length(truebits)
    output = output + truebits(m)*str2num(dn(end-m+1));
end

end
%--------------------------------------------------------------------------