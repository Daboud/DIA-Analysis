%Final data is output as vector "clipour", and automatical copied to user's
%clipboard

klear=1; % set to 0 to keep previous frame measurements and change only the
         % final analysis
if klear==1;
    clear
   klear=1;
end

%Setting anumber of useful parameter that vary by video
vt=1.74; %tangential velocity of the drop, from another script
calib=0.0334% calibration const. of video, mm/pix
adjustfac=1.0; % adjustment factor, separates the measurement line from the 
               % regressed position line of the surface by a pixel or two
vs=0; %surface velocity, from another script
frate=20; %framerate of the video (kiloframes/s)
bal=20; % Threshold brightness value to separate drop and surface from back
        % ground
fskip=10; %Show user an image of every Xth frame and measurement
ignoreYs=[] % Y-positions in the video to ignore in measurement
ignoreIs=[] %frames to ignore due to video errors
format short

%turns off a repetative MATLAB error message
id =['MATLAB:colon:operandsNotRealScalar']; 
warning('off',id);

%Obtaining the total number of frames of the video
files = dir('*.tif');
file=files(end,1); %all files in folder
file=file.name; %name of file item in folder
f2=str2num(file(7:9)); %ending frame
f1=1; %start from frame 1

if klear==1
xy2=[]; %empty matrices for later
xy1=[];

colormap(gray); %get figure position and colour scheme prepared
MainFig = figure;
  set(MainFig, 'Position', [100, 200, 1120, 430]);
close(1)

for i=1:f2; %begin iteration of frame measurements
disp(['Frame ',num2str(i)])
if i<10 %call on correct image file for current frame
    I0 = double(imread(['vid00000',num2str(i),'.tif']))/4;
elseif i<100
    I0 = double(imread(['vid0000',num2str(i),'.tif']))/4;
elseif i>99
    I0 = double(imread(['vid000',num2str(i),'.tif']))/4;
end

I0=I0(:,:,1); %remove double layers. I0 is original unprocessed image
[m,n]=size(I0); %height, width of images

%turn black & white
I=I0; %I will be black/white version of image
black=find(I<bal);
I(black)=0;
white=find(I>=bal);
I(white)=256; 
I=256-I;
I=imfill(I);   %Image filling, remove white glare inside droplet
I=256-I;

%trunkL = "left truncate", essentially test points further inward from the 
% boundaries of the image, usually because the edge of the sample surface 
% comes into view
trunkL=1;
if i>40 %in case surface pops into view at certain frame
    trunkL=round(((i-40)/frate)*vs/calib)+1; %X axis line of truncation
end
trunkR=n;
% trunkR=150+round(((i)/frate)*vs/calib);
if trunkR>n
    trunkR=n;
end
% trunkL=1;
% trunkR=n;

% find the border between surface and background at each side of the image
% to estimate surface position line (green)
xx=[]; %collect points at each end
yy=[];
for j=1:6
    Y1=j; %for the top 6 rows
    X1=min(find(I(Y1,trunkL:trunkR)==256))-1+trunkL; %find the lowest x value that is white
if X1==trunkL %in case surface hits left wall
    X1=j+trunkL; %for the leftmost six columns
    Y1=length(find(I(:,X1)==256))+1; %find the length of y locations that are white
end
Y2=m-j+1; %for the bottom six rows
X2=find(diff(I(Y2,1:trunkR)==0)==-1); %find where the background starts
if length(X2)==0%length(I(Y2,1:trunkR)) %in case line hits the right wall
    X2=trunkR-j+1; %for the rightmost six columns
    Y2=m-min(find(flipud(I(:,X2))==256))+2; %find height of surface
end
xx=[xx;X1;X2]; %accumulated vector of surface positions
yy=[yy;Y1;Y2];
end
slashco=polyfit(xx,yy,1); %estimate of surface position (slash coefficient)
ms=slashco(1); %slope angle in ratio dy/dx
ang=atan(ms); %angle in radians
coefsG=slashco-[0 adjustfac/cos(ang)]; %coefs of the measurement line,
% above which all black pixels are assumed to be droplet, not surface.
%coefsG = "coefficients of green line"
xline=[0 n]; %for plotting measurement line later
yline=polyval(coefsG,xline);

%Im = "image Minus surface", essentially removing all traces of surface
%and leaving only the drop in the image for analysis
Im=I;
for x=1:n %for each x location
y=coefsG(1)*x+coefsG(2); %find the surface's y location
y=round(y);
if y>2 
    Im(round(y):end,x)=256; %and turn it all white below
else
    Im(:,x)=256; %%left hand side bar, just in case (avoiding error)
end
end
Im(ignoreYs,:)=256; %Ignore any points within a certain range

[Y1 X1]=find(Im==0); %find where the drop is
X1=mean(X1);
Y1=mean(Y1); %get avg x-y positions
xg=(Y1+1/ms*X1-coefsG(2))/(ms+1/ms); %project onto green line, essentially
%keeping things 1-dimensional
yg=polyval(coefsG,xg); % now we have x-y location of the projected drop

%Dynamic tracking throughout impact process, over all i's
xy2=[xy2; xg yg]; %xy2 tracks position of drop from 1st to last frame
%for analysis of its dynamic movement (red circle in image)
if i==1;
    yg1=yg; %on first frame, remember the y-position of impact
end
%xy1 tracks the initial impact center as it moves with the surface
y1=yg1; %its y-position never changes
x1=(y1-coefsG(2))/ms; %but the x-position slides with the surface
xy1=[xy1;x1 y1]; %tracking impact center, the blue square

%every X frames, pause and show a plot. X defined by fskip
if i==1 | round(i/fskip)==i/fskip 
subplot(1,2,1)
image(I0) %original image on left
hold on
plot(xline,yline,'--g') %measurement line
plot(xx,yy,'om') %surface line regression points
subplot(1,2,2)
image(Im) %adapted image on right
colormap(gray)
hold on
plot(xline,yline,'--g')%measurement line
plot(xg,yg,'ro','LineWidth',1.5,'MarkerSize',7) %projected drop footprint
plot(x1,y1,'bs','LineWidth',1.5,'MarkerSize',7) %impact center
pause
end
end
end %here ends the iteration, and begins final analysis

close %previous plot

d=xy1-xy2; % pixel distance between red circle and blue squre, based on 
           % green line projection, over all i's
dist=sqrt(d(:,1).^2+d(:,2).^2)*calib; %dynamic mm distance
t=[[1:length(dist)]/20-.05]'; %time vector
t(ignoreIs)=[]; %in case of bad or glitchy frames
dist(ignoreIs)=[];
L0=(vt*t/1000+4.9*sin(ang)*(t/1000).^2)*1000; %ideal sliding length based
% on tangential velocity and gravitation effect

x=t; %easier variables to use for regression
yr=dist;
% regression, finding 95% CI of L/L0
[b, bint] = regress(yr(:), [x(:) ones(numel(x),1)]);

%Get the name of the original file to label graph
CurDir=pwd;
Slass=find(CurDir=='\');
VName=CurDir(max(Slass)+1:end);
% title(VName);

%final plot of dynamic drop movement across surfaces vs. L0
plot(t,L0,'--','LineWidth',2,'Color',[0 0.7 0])
hold on
plot(t,dist,'k.','MarkerSize',9)
plot(t,polyval(b,t),'r','LineWidth',1)
xlabel ('Time (ms)')
ylabel ('Length (mm)')
title(VName)
legend('L_0','L (measured)','L (fit)','Location','SouthEast')
set(gca, 'box', 'off')
set(gca,'fontsize',20)

%get final values of L/L0 and the 95% CI, and copy to clipboard
L=b(1)*t(end);
dL=(bint(3)-bint(1))*t(end);
clipout=[L,dL]
clipout2=num2str(clipout);
clipboard('copy',clipout2)

%print the final plot
name=['LL0 Dynamic'];
prt=1; %switch this function on/off
if prt==1
width = 7.2;     % Width in inches
height = 5.4;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(name,'-dpng','-r300');% Save the file as PNG
end