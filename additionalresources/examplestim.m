function examplestim

clear
close all;

E.levelsC5 = [48 0 48 48]./100;
E.levelsC7 = [0 48 48 48]./100;

stimorient = 360.*(rand(36,1));

ST.npixelsperdegree = 36;       % VPixx pixel resolution at 57cm

ST.gratingsize = ST.npixelsperdegree*2;         %
ST.gratingradiusfactor = ST.npixelsperdegree*2;  % previously 128
ST.SF = 2;
ST.ncycles = ST.SF*ST.gratingsize/ST.npixelsperdegree;

screens=Screen('Screens');
screenNumber=max(screens);
rect = [1 1 1024 1024];
[w, winRect] = Screen('OpenWindow',screenNumber,128,rect);
ST.greylevel = 128;
[width, height] = Screen('WindowSize', w);

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

ST.dotlevels = round(255*rand(3,ndots));

ST.dotsize = 8;
ST.centre = [width/2, height/2];

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

comp = sgrating.*window.*E.levelsC5(1);
comp = (1+comp)/2;
targettexture = Screen('MakeTexture', w, comp, [], [], 2);


comp = sgrating.*window.*E.levelsC5(1) + rot90(sgrating).*window.*E.levelsC7(2);
comp = (1+comp)/2;
plaidtexture = Screen('MakeTexture', w, comp, [], [], 2);


Screen('FillRect',w, 128);

Screen('DrawTextures', w, targettexture, [], rectlist, stimorient+45);

Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots
Screen('Flip', w);

im = Screen('GetImage',w);

WaitSecs(2);

im = im(256:1791,256:1791,:);
imwrite(im, 'gratinggrid.jpg','jpeg');


Screen('FillRect',w, 128);

Screen('DrawTextures', w, targettexture, [], rectlist, stimorient-45);

Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots
Screen('Flip', w);

im = Screen('GetImage',w);

WaitSecs(2);

im = im(256:1791,256:1791,:);
imwrite(im, 'gratinggrid2.jpg','jpeg');



Screen('FillRect',w, 128);


Screen('DrawTextures', w, plaidtexture, [], rectlist, stimorient+45);

Screen('DrawDots', w, ST.dotcoords, ST.dotsize, ST.dotlevels, ST.centre, 0);  % change 1 to 0 to draw square dots
Screen('Flip', w);

im = Screen('GetImage',w);

WaitSecs(2);

im = im(256:1791,256:1791,:);
imwrite(im, 'plaidgrid.jpg','jpeg');



Screen('CloseAll');




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