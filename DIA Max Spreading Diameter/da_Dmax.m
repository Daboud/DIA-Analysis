%Final data is output as an array called "clipout"

klear=1; % set to 0 to keep previous frame measurements and change only the
         % final analysis
if klear==1;
    clear
   klear=1;
end

%-------Setting anumber of useful parameter that vary by video----------
format compact
calib=0.0334;% calibration const. of video, mm/pix
adjustfac=2.0; %separation between surface line and measurement line
adjustfacM=2.0 ; %same as above, but specific to the maximum diameter
swagfac=0.3; %this 'swaggers' the measurement line a little, to detect
%the most inward position of the 3-phase contact line
vs=2; %velocity of the surface, from other code
bal=15; % Threshold brightness value to separate drop and surface from back
        % ground, for contact diameter (Dc) measurement
bal2=15; %same, but for maximum diameter (Dm) measurement
fskip=3; %show user images of the analysis every X frames
ignoreYs=[]; %remove a section of the image from analysis
frate=20; % framerate of video in kiloframes/s
twist=1+0/100; %rotates the surface line slightly in case of distortion

%------determines the total frame count from the root folder--------
f1=1;  
files = dir('*.tif');
file=files(end,1); %all files in folder
file=file.name; %name of file item in folder
f2=str2num(file(7:9)); %ending frame

id =['MATLAB:polyfit:RepeatedPointsOrRescale']; %remove a pesky error msg
warning('off',id);

%--------------------Should the function...-------------------------
prt=0; %print images of each frame to create a marked video later (SI in
%2018 maximum diameter and sliding length paper)
plotMuphill=1;%plot max diameter wrt impact center?
plotCuphill=1;%plot contact diameter wrt impact center?
paus=1; %pause periodically to view images?
plotDvsT=1; %plot contact diameter overall?

colormap(gray);
dmp=[]; %store progressive Dm
dcp=[]; %store progressive Dc

FigHandle = figure;%set figure size and position for later
  set(FigHandle, 'Position', [50, 50, 1300, 630]);
close(1)

%----------------Begin iteration of measurements for each frame-----------
clipout=[calib bal adjustfac vs 0 0 0];
if klear==1
for i=f1:f2
disp(['Frame ',num2str(i)])
%import iage for frame i
if i<10
    I0 = double(imread(['vid00000',num2str(i),'.tif']))/4;
elseif i<100
    I0 = double(imread(['vid0000',num2str(i),'.tif']))/4;
elseif i>99
    I0 = double(imread(['vid000',num2str(i),'.tif']))/4;
end
I0=I0(:,:,1);
[m,n]=size(I0); %height, width of image

%useful for adjusting analysis parameters as surface slides across screen
%and may be distorted
if i>40
   %bal=7
   %fskip=1;
   %adjustfac=1.6;
   %swagfac=0.4;
end

%-------black & white version of original image, from threshold (bal)------
I=I0;
black=find(I<bal);
I(black)=0;
white=find(I>=bal);
I(white)=256; %turn black & white

%same as above, for Dm
I2=I0;
black=find(I2<bal2);
I2(black)=0;
white=find(I2>=bal2);
I2(white)=256; %turn black & white

%----------determination of surface line, linear regression--------------
%trunkL = "left truncate", essentially test points further inward from the 
% boundaries of the image, usually because the edge of the sample surface 
% comes into view
trunkL=1; % 
if i>48 %in case surface left edge pops into view at certain frame
    trunkL=round(((i-47)/frate)*vs/calib)+1; %X axis line of truncation
end
trunkR=140+round(((i)/frate)*vs/calib);
% trunkR=n;
if trunkR>n;trunkR=n;end

% find the border between surface and background at each side of the image
% to estimate surface position line (green)
xx=[]; %collection of points found
yy=[];
for j=1:6
    Y1=j; %for the top 6 rows
    X1=min(find(I(Y1,trunkL:trunkR)==256))-1+trunkL; %find the lowest white 
    %x value
if X1==trunkL %in case surface hits left wall
    X1=j+trunkL; %for the leftmost six columns
    Y1=length(find(I(:,X1)==256))+1; %find the length of white y locations
end
Y2=m-j+1; %for the bottom six rows
X2=min(find(diff(I(Y2,1:trunkR)==0)==-1));%how long do the black pixels go?
if length(X2)==0%length(I(Y2,1:trunkR)) %in case line hits the right wall
    X2=trunkR-j+1; %for the rightmost six columns
    Y2=m-min(find(flipud(I(:,X2))==256))+2;
end
xx=[xx;X1;X2]; %accumulate points found for j=1:6
yy=[yy;Y1;Y2];

end
slashco=polyfit(xx,yy,1); %"slash coefficient"=surface line
slashco(1)=slashco(1)*twist;%twisted slightly if neccessary
ms=slashco(1); %slope angle in ratio dy/dx
ang=atan(ms); % slope in radians
coefsG=slashco-[0 adjustfac/cos(ang)]; %coefs of contact line for Dc
coefsM=slashco-[0 adjustfacM/cos(ang)]; %coefs of line for Dm

%-----------Remove surface from "Im" for measurement of Dm------------
Im=I2; %Matrix to find max diameter points
for x=1:n %for each x location
y=coefsM(1)*x+coefsM(2); %find the slope's y location
y=round(y);
if y>2
    Im(round(y):end,x)=256; %and turn it all white below
else
    Im(:,x)=256; %left hand side bar, just in case
end
end
Im(ignoreYs,:)=256; %Ignore any points within a certain range

%-----------plot the original image, the Dm image, and the Dc image
subplot(2,3,1)
colormap(gray);
image(I0) %original
hold on
plot(polyval(coefsG,1:n),'--g');%measurement line
plot(polyval(slashco,1:n),'--r') %surface line
subplot(2,3,[2 3 5 6])
image(I) %Dc image
hold on
plot(polyval(coefsG,1:n),'--g');
plot(polyval(slashco,1:n),'--r');
plot(xx,yy,'or') %surface tracking markers
subplot(2,3,4)
image(Im) %Dm image
hold on



%-----------------------------Dc measurement----------------------
%This code 'swaggers' the measurment line slightly (swagfac) to see if a
%more inwards point can be found just above/below the measurement line.
%This helps when the focus or lighting affects one side of the drop's
%contact line
leftbasket=[];%"basket" is a collection of measured points. The code will
%later select the most inwards Left side.
rightbasket=[]; % right side
for swag=0:4; % a total of 5 times
    mtp=swag-2; % move meas line slightly above or below original location
    %coefs of adjusted line:
    coefs2=slashco-[0 adjustfac/cos(ang)]-mtp*[0 swagfac/cos(ang)]; 
    if swag==0 | swag==4% plot in yellow furthest extends of the swagger
    subplot(2,3,[2 3 5 6]);
    plot(polyval(coefs2,1:n),'--y')
    end
    
line=[];%through swaggered measurement line, record brightness
pix=1; %counter for iterations below
for jj=1:cos(ang):n 
    x=round(jj);%for each x location
    y=coefs2(1)*jj+coefs2(2); %follow measurement line
    y=round(y); % y value at x
    if y<1 %if y is outside image
        line(pix)=256; % background
    elseif y>m %if y is outside image
        line(pix)=256; %background
    elseif y<=m %if y is within image
    line(pix)=I(y,x); %record the brighness value
        if sum(y==ignoreYs)>0  %if the tested Y loc is in ignore range
            line(pix)=256; %then just make it white on the line
        end
    end
    pix=pix+1;
end
lmin=min(find(line<bal)); %leftmost footprint of drop
lmax=max(find(line<bal))+1; % rightmost footprint of drop
%x-y locations of drop's leftmost footprint (1) and rightmost (2)
x1=lmin*cos(ang);%left X actual pixel locations
y1=coefs2(1)*lmin*cos(ang)+coefs2(2); 
x2=lmax*cos(ang); %right X actual pixel locations
y2=coefs2(1)*lmax*cos(ang)+coefs2(2);
%collect all 5 swaggers, on left and right
leftbasket=[leftbasket; x1 y1];
rightbasket=[rightbasket; x2 y2];
end

%Finding which of the 5 collected points is most inwards
m2=-1/ms; %inverted slope, perpendicular to measurement line
b1=coefsG(2); %y=mx+b coefficient of measurement line
m1=coefsG(1);

for g=1:5; %for each of five swaggered points on the left
x2=leftbasket(g,1); %get the x-y coordinates
y2=leftbasket(g,2); 
b2=y2-m2*x2;%project it onto the measurement line
x=(b1-b2)/(m2-m1);
lefty(g)=m1*x+b1;  %and get the y-intercept
end

for g=1:5; %same idea for right side of the drop
x2=rightbasket(g,1);
y2=rightbasket(g,2); 
b2=y2-m2*x2;
x=(b1-b2)/(m2-m1);
righty(g)=m1*x+b1;
end

n1=find(lefty==max(lefty));n1=n1(1); %of the 5 points collected, which one
% has the most inwards loction projected onto the measurement line?
n2=find(righty==min(righty));n2=n2(1);

x01=leftbasket(n1,1); %decided most inwards, considering all 5 swaggered
%points
y01=leftbasket(n1,2);
x02=rightbasket(n2,1);
y02=rightbasket(n2,2);

x1F=(y01+1/ms*x01-coefsG(2))/(ms+1/ms); %projecting them onto meas line
y1F=polyval(coefsG,x1F);
x2F=(y02+1/ms*x02-coefsG(2))/(ms+1/ms);
y2F=polyval(coefsG,x2F);

%If the sample is flat, none of this works well. So this is a quick fix
if abs(ms)<0.1
    x1F=max(leftbasket(:,1));
    y1F=polyval(coefsG,x1F);
    x2F=min(rightbasket(:,1));
    y2F=polyval(coefsG,x2F);
end

%In the regular pop-up graph to track progress, show red Xs for the chosen
%swaggered points, and magenta Xs for their projected locations on the
%measurement line
subplot(2,3,[2 3 5 6]);
plot(x01, y01,'rx','LineWidth',1)
plot(x02, y02,'rx','LineWidth',1)
plot(x1F, y1F,'mx','LineWidth',2,'MarkerSize',8)
plot(x2F, y2F,'mx','LineWidth',2,'MarkerSize',8)
%------Here ends the lengthy 'swagger' process---------------

%on the first frame, define the impact CENTer between L/R footprints
if i==1;
    ycent=(y1F+y2F)/2;
    xcent=(x1F+x2F)/2;
end
%Henceforth, X-center has to slide along with the surface
if vs>0
    xcent=(ycent-coefsG(2))/coefsG(1);
end
subplot(2,3,[2 3 5 6]);plot(xcent,ycent,'bs','MarkerSize',8,...
    'LineWidth',1.5);

%This measures the contact positions of the front and tail of the lamella
%with respect to the impact center. Eventually used in the 'uphill' plots
XfrontCp(i)=(y1F-ycent)/sin(ang); %length in pixels
XtailCp(i)=(y2F-ycent)/sin(ang);
if ang<pi/4 %more accurate for flatter surface tilt angles
    XfrontCp(i)=(x1F-xcent)/cos(ang);
    XtailCp(i)=(x2F-xcent)/cos(ang);
end

dcp(i)=-XfrontCp(i)+XtailCp(i); % Dc of the drop's footprint in pixels
%-----------------------END CONTACT DIAMETER measurement------------------

%----------------------START MAX DIAMETER MEASUREMENT---------------------
[Imby Imbx]=find (Im==0); %find all black pixels (drop) in the image

m2=-1/ms; %inverted slope, perpendicular to measurement line
b1=coefsG(2); %y=mx+b coefficient of measurement line
m1=coefsG(1);

ygrn=[];%??
for g=1:length(Imbx); %for each pixel
x2=Imbx(g); % get the x-y coordinates
y2=Imby(g); 
b2=y2-m2*x2; %project them onto the measurement line
x=(b1-b2)/(m2-m1); %and get the prjected x-coordinate value
ygrn(g)=m1*x+b1; %and projcted y-coordinate value
end

n1=find(ygrn==min(ygrn));%select the most inwards and outwards points by 
%finding which points in the ygrn vector have the highest/lowest values
n2=find(ygrn==max(ygrn));
n1=n1(1);
n2=n2(1);

subplot(2,3,4) %mark the corresponding pixels with red circles in the image
plot(Imbx(n1),Imby(n1),'or')
plot(Imbx(n2),Imby(n2),'or')

y1=min(ygrn);%y-coordinates of selcted points
y2=max(ygrn);

XfrontMp(i)=(y1-ycent)/sin(ang); %x-coordinates of selcted points, wrt the
%impact center
XtailMp(i)=(y2-ycent)/sin(ang);

if ang<0.1 %for very low slope, abandon fancy analysis and pick outmost pix
    XfrontMp(i)=min(Imbx)-xcent;
    XtailMp(i)=max(Imbx)-xcent;
end

dmp(i)=-XfrontMp(i)+XtailMp(i); %measured Dm value in pixels
%------------Here ends the maximum diameter measurement-----------------

%This is useful for making marked frames to build a cool video
if prt==1
pf = figure;
  set(pf, 'Position', [50, 50, 1300, 900]); 
image(I0) 
hold on
colormap(gray);
p(1)=plot(x1F, y1F,'ro','LineWidth',1,'MarkerSize',8);
p(2)=plot(x2F, y2F,'ro','LineWidth',1,'MarkerSize',8);
p(3)=plot(xcent,ycent,'s','LineWidth',1,'MarkerSize',9,'Color',[0 0.5 1]);
axis off
print(num2str(i),'-dpng','-r300');
close
end

%pause periodically (every fskip frames) the show user progress
if ((i)/fskip==round(i/fskip) | i==1) & paus==1; 
pause
elseif i==f2
    pause
end

end
end

dc=dcp*calib; %Dc in mm instead of pixels
dm=dmp*calib; %Dm
t=(f1:f2)/frate-1/frate;%time in ms starting from 0
%accumulatin of measured data
data=[t' XfrontCp'*calib XtailCp'*calib dc' XfrontMp'*calib ...
    XtailMp'*calib dm'];
clipout=[clipout; data]; % cap that with chosen settings for repeatability
close

if plotDvsT==1 %plot of Dc and Dm versus time
figure
plot(t,[dc',dm'],'x','MarkerSize',4) 
end

if plotMuphill==1 %maximum extension tracking wrt impact center
figure
plot(t,XfrontMp*calib,'or','MarkerSize',4) 
hold on
plot(t,XtailMp*calib,'or','MarkerSize',4,'MarkerFaceColor','r') 
end

if plotCuphill==1 % 3-phase contact line tracking wrt impact center
figure
plot(t,XfrontCp*calib,'or','MarkerSize',4)
hold on
plot(t,XtailCp*calib,'or','MarkerSize',4,'MarkerFaceColor','r') 
end

disp(['Ignoring Y from:  ',num2str(min(ignoreYs)),' to '...
    num2str(max(ignoreYs))])

