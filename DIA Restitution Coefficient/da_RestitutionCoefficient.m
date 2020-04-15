%The output of this program will be in a vector called "clipout"

klear=1;% clear old data and re-process? Turn off if you need to adjust
%prameters like vs and theta but don't want to reprocess frames
if klear==1;
    clear
    klear=1;
end

%Some important switches and parameters for analysis
regress=1; %regress and find RC?
prt=0;%Print frames to make figures?
imp=1; %Should the dY of the drop be tracked by this code?
dYvid=0;%Or should it be entered manually by the user?
adjustfac=3;%Used for impact center location finding
theta=25.28; %degree of surf tilt
calib=0.0334;% calibration const of video, in mm/pix
frameskip=10;%visualize progress every X frames
kills=[];%kill data at anomalous data points, if code misinterprets video
vs=2.175101194; %surface velocity
bal=15;%Threshold brightness, to separate drop and surface from background
fr=20; %framerate of video, kiloframes/s
afterlimit=170; %maximum number of frames after drop detachment to analyze


id =['MATLAB:colon:operandsNotRealScalar']; %turn off repetetive error msg
warning('off',id);

%This snippet finds how many frames are in the video in the working folder
files = dir('*.tif');
file=files(end,1); %all files in folder
file=file.name; %name of file item in folder
f2=str2num(file(7:9)); %ending frame
f1=1; %start from frame 1

if klear==1; %is droplet tracking performed?
swch=0; %This 'switch' helps later, determining impact/detachment moments
Isend=[];%This helps the code keep record of the previous frame, to determ-
%ine the moment of detachment

colormap(gray);%Setting position and height of figure for later
MainFig = figure;
  set(MainFig, 'Position', [250, 50, 1120, 330]);
close(1)

for i=(f1):f2 %from first to last frame, track drop location
disp(['Frame ',num2str(i)]) %for user
if i<10 %calling the correct image file for frame i. I0=original image
    I0 = double(imread(['vid00000',num2str(i),'.tif']))/4;
elseif i<100
    I0 = double(imread(['vid0000',num2str(i),'.tif']))/4;
elseif i>99
    I0 = double(imread(['vid000',num2str(i),'.tif']))/4;
end
[m,n]=size(I0); %height, width of image

I=I0;%I=black&white version of image
black=find(I<bal);
I(black)=0;
white=find(I>=bal);
I(white)=256; 
I=256-I;
I=imfill(I); %Image filling, removes white glare spots inside droplet
I=256-I;
I(1:m,end)=0; I(end,1:n)=0;%This adds a black border to the bottom and 
%right of the adapted image. It does not interefere with analysis, but
%helps the code ignore any lens effects in the corners of the video

%This snippet identifies black sections in the image, and keeps only items
%within a certain size (total pixel count). Removes surface, keeps drop!
Ibin=(256-I)/256; %convert I to binary information
LB = 100; %smallest blob to keep
UB = 3000; %largest blob allowable
Iout = xor(bwareaopen(Ibin,LB),  bwareaopen(Ibin,UB)); %keep only between
I2=(1-Iout)*256; %create image using only kept (drops)


%--------------------------WEIGHTED DROPS--------------------------
%This is needed when the rebounded drop fragment into pieces. This code
%identifies each piece, estimates their size based on their 2D silhouette,
%and then locates the center of gravity to track the momentum of the liquid
Iw=1-I2/256;%image of just the drop, in Binary form
CC = bwconncomp(Iw); %CC has a bunch of wrapped up info about sizes
numPix = cellfun(@numel,CC.PixelIdxList); %Number of pixels of each item
numObject = CC.NumObjects; %total number of items

if numObject==0; %numObject=0 implies surf-liq contact.
    disp('Contact') %info for user
    xmp(i)=mean([]);%sets drop x-loaction to NaN while in contact
    ymp(i)=mean([]); 
else
    
exes=[]; %used below
for j=1:numObject %for each drop fragment
Ij=ones(m,n); %create white image
Ij(CC.PixelIdxList{j})=0; %set drop area black
[y,x]=find(Ij==0); %find black pixels
xtemp=mean(x); %get mean coordinates
ytemp=mean(y);
exes=[exes;xtemp ytemp numPix(j)^1.5]; %this matrix funnels the data down
%to the x and y mean location for each drop, and its weigh on the center-
%of gravity calculation, based on its crossectional area
if numPix(j)==(m+n-1)
    exes(end)=0 %ignores the black border around the image!
end
end

%This finds the center of gravity
SumWeight=sum(exes(:,3)); %full weight
exes0=exes;%reserve copy for later
exes(:,3)=exes(:,3)/SumWeight; %turns right column to weight fractions
exes(:,1)=exes(:,1).*exes(:,3);%weight affects X and Y columns (1,2)
exes(:,2)=exes(:,2).*exes(:,3);
xmp(i)=sum(exes(:,1)); %sum column to get weighted position
ymp(i)=sum(exes(:,2));
end
%--------------------------END WEIGHTED DROPS-------------------------

%Every so many frames, pause and plot progress for user
if i==1 | round(i/frameskip)==i/frameskip | i==f2
subplot(1,3,1)
image(I0)
hold on
subplot(1,3,2)
image(I)
hold on
subplot(1,3,3)
image(I2)
colormap(gray)
hold on
for j=1:numObject
plot(exes0(j,1),exes0(j,2),'or','LineWidth',1)
plot([exes0(j,1),xmp(end)],[exes0(j,2),ymp(end)],'g-','LineWidth',1)
end
plot(xmp(end),ymp(end),'+b','MarkerSize',20,'LineWidth',1)
pause
end

%------------START IMPACT CENTER DETECTION TO DETERMINE dY-------------
%Because gravity affects the restitution coef, this section finds the
%drop's impact center on the surface and the leave-off point, so that the
%change in height can be measured, and the change in gravitational
%potential energy can be used to adjust the calculated RC accordingly. The
%measured change in height is dY
if imp==1 %turn off if using manually entered dY
IsendPr=Isend; %save copy of previous image to find detachment time
Isend=I(1:m-1,1:n-1); %Remove black border
[m2,n2]=size(Isend); %height, width of image used for impact center finding
if (xmp(i)>0)==0 & swch==0 %This condition detects onset of contact
swch=swch+1; %trigger switch, indicating contact started

%This snippet estimates the surface position using found points and regress
xx=[]; %found points at each end
yy=[];
for j=1:6
    y1=j; %left side
    x1=min(find(Isend(y1,:)==256))-1;
if x1==0 %in case angle hits left wall
    x1=j;
    y1=length(find(Isend(:,x1)==256))+1;
end
y2=m2-j+1; %right side
x2=length(find(Isend(y2,:)==0)); %how many black pixels
if x2==n2 %in case line hits bottom wall
    x2=n2-j+1;
    y2=length(find(Isend(:,x2)==256));
    y2=m2-min(find(flipud(Isend(:,x2))==256))+1;
end
xx=[xx;x1;x2]; %accumulate points j times for 2j long vector
yy=[yy;y1;y2];
end
coL=polyfit(xx,yy,1); %coefficients of surface line
xl=[1 n];%Plot coL on graph for user
yl=xl.*coL(1)+coL(2);
ms=coL(1); %slope angle in ratio dy/dx
ang=atan(ms); %state in radians

%Get coefs of measurement line, adjusted from surface line by adjustfac
coefs2=coL-[0 adjustfac/cos(ang)]; %shift the slash line x pixels away from surface

%along the line of contact between surface and drop, record the brightness
%value, select the black pixels (drop loacation), and estimate impact centr
line=[];
pix=1;
for j=1:cos(ang):n-1 %for each x location, but adjusted to not skip too 
    %fast on high angles!!
    x=round(j); %discrete x value
    y=coefs2(1)*j+coefs2(2); %follow line through bottom of drop
    y=round(y); %discrete y value
    if y<1 %avoid errors
        line(pix)=256;
    elseif y>m-1 % avoid errors
        line(pix)=256;
    elseif y<=m-1
    line(pix)=Isend(y,x); %and make a vector of all the shades it crosses
    end
    pix=pix+1;
end
ys=coefs2(1)*find(line==0)*cos(ang)+coefs2(2); %y locations of black pixels
xs=(ys-coefs2(2))/coefs2(1);%corresponding x-locations
ximp=mean(xs); %impact center locations
yimp=mean(ys);

%To check if the impact center was estimated correctly
figure
colormap('gray')
image(Isend)
hold on
plot(xl,yl,'--r')%surface line
plot(xx,yy,'or') %used to defined surface line
plot(ximp,yimp,'go') %impact center
plot(polyval(coefs2,1:n),'--b');% measurement line
pause
close

Ygrav1=yimp  %y-position at impact
end

if (xmp(i)>0)==1 & swch>0
    swch=-1;

%SAME AS ABOVE BLOCK OF CODE-----------------------------------------------
%This snippet estimates the surface position using found points and regress
xx=[]; %found points at each end
yy=[];
for j=1:6
    y1=j; %left side
    x1=min(find(Isend(y1,:)==256))-1;
if x1==0 %in case angle hits left wall
    x1=j;
    y1=length(find(Isend(:,x1)==256))+1;
end
y2=m2-j+1; %right side
x2=length(find(Isend(y2,:)==0)); %how many black pixels
if x2==n2 %in case line hits bottom wall
    x2=n2-j+1;
    y2=length(find(Isend(:,x2)==256));
    y2=m2-min(find(flipud(Isend(:,x2))==256))+1;
end
xx=[xx;x1;x2]; %accumulate points j times for 2j long vector
yy=[yy;y1;y2];
end
coL=polyfit(xx,yy,1); %coefficients of surface line
xl=[1 n];%Plot coL on graph for user
yl=xl.*coL(1)+coL(2);
ms=coL(1); %slope angle in ratio dy/dx
ang=atan(ms); %state in radians

%Get coefs of measurement line, adjusted from surface line by adjustfac
coefs2=coL-[0 adjustfac/cos(ang)]; %shift the slash line x pixels away from surface

%along the line of contact between surface and drop, record the brightness
%value, select the black pixels (drop loacation), and estimate impact centr
line=[];
pix=1;
for j=1:cos(ang):n-1 %for each x location, but adjusted to not skip too 
    %fast on high angles!!
    x=round(j); %discrete x value
    y=coefs2(1)*j+coefs2(2); %follow line through bottom of drop
    y=round(y); %discrete y value
    if y<1 %avoid errors
        line(pix)=256;
    elseif y>m-1 % avoid errors
        line(pix)=256;
    elseif y<=m-1
    line(pix)=Isend(y,x); %and make a vector of all the shades it crosses
    end
    pix=pix+1;
end
ys=coefs2(1)*find(line==0)*cos(ang)+coefs2(2); %y locations of black pixels
xs=(ys-coefs2(2))/coefs2(1);%corresponding x-locations
ximp=mean(xs); %impact center locations
yimp=mean(ys);

%To check if the impact center was estimated correctly
figure
colormap('gray')
image(Isend)
hold on
plot(xl,yl,'--r')%surface line
plot(xx,yy,'or') %used to defined surface line
plot(ximp,yimp,'go') %impact center
plot(polyval(coefs2,1:n),'--b');% measurement line
pause%--------------------------------------------------------------------
close

    Ygrav2=yimp %y-position after detach
    pause
    dym=(Ygrav2-Ygrav1)*calib; %final dY value
end

else
dym=dYvid*calib; %If manually entering dY instead of detection
end
%------------END IMPACT CENTER DETECTION TO DETERMINE dY-------------

%This section prints images of the analysis to make figures
if prt==1
    if round((i-1)/15)==((i-1)/15)
pf = figure;
  set(pf, 'Position', [50, 50, 1300, 300]);
  figure(1)
image(I) 
hold on
colormap(gray);
plot(xmp(end),ymp(end),'xr','MarkerSize',8,'LineWidth',1)
set(gca,'position',[0 0 1 1],'units','normalized')
axis off
% p(2).LineWidth = 10
print([num2str(i)],'-dpng','-r300');
close
    end
end

end %End frame-by-frame analysis
close%various figures
close
close

xm=xmp*calib;%x-position vector in mm 
ym=ymp*calib;%y-position vector in mm

end %droplet tracking performed?

if regress==1; %regress and find RC?
    
%Just below finds the impact and detach frames of drop contact
a=(xm>0); %where did the tracking fn not just error?
b=diff(a); 
c=find(b~=0); %find spots where flips from zero fields / NaN to numbers
stnan=min(c); %touch surface frame
endnan=max(c)+1; %leave surface frame

%------------Linear regression to get drop velocities----------------
%After Impact...
xmaft=xm(endnan:end); %take all spots after detach
ymaft=m*calib-ym(endnan:end); %flip it
if (length(xm)-endnan)>afterlimit%analyze only for so many frames
    xmaft=xmaft(1:afterlimit); 
    ymaft=ymaft(1:afterlimit);
end
fl2=length(xmaft); %length of the time period after detach in frames
t2=[1:fl2]/fr; %create time vector for detachment tracking
xmaft(kills)=[];%Remove data in glitchy spots
ymaft(kills)=[];
t2(kills)=[];

%This is a manual way to essentially do 2nd order curve fitting with the
%second order term set at -9.8 for gravity. So this snippet alters the
%tracked y-position of the drop to what it would have been without gravity
%interfering
dvy=t2*-9.8;  %Account for change in speed from gravity
ymaftadj=ymaft-.5*dvy.*t2/1000;  %Account for loss of distance by gravity
yfit2=polyfit(t2,ymaftadj,1); %linear fit for altered y-position
y2 = polyval(yfit2,t2); %fitting curve defined
xfit2=polyfit(t2,xmaft,1); %linear fit for x
x2 = polyval(xfit2,t2); %fitting curve defined

%Before impact...
xmbef=xm(1:stnan); %take all spots before contact
ymbef=m*calib-ym(1:stnan); %flip it
fl1=length(xmbef); %length of the time period
t1=[-fl1+1:0]/fr; %start from t=0 for this case
yfit1=polyfit(t1,ymbef,2); %fit quadratic
y1 = polyval(yfit1,t1); %fitting curve defined for drop falling onto surf

%---------Plotting all three measured velocities---------------------
FigHandle = figure; %Open figure window to right size
  set(FigHandle, 'Position', [10, 50, 1330, 500]);%(959-567)/n*m]);

subplot(1,3,1) %vy1 before touch
plot(t1,y1,'LineWidth',6,'Color',[0 .8  .8])
title(['Impacting y-velocity';
       '                    '])
xlabel('Time (ms)')
ylabel('Position (mm)')
hold on
plot(t1,ymbef,'ow','MarkerSize',4,'MarkerFaceColor','k','LineWidth',0.25)
legend('2^n^d Order Fit','Data','Location','SouthWest')
set(gca, 'box', 'off')
set(gca,'fontsize',20)
hold off

subplot(1,3,2) %vx after detach
plot(t2,x2,'LineWidth',6,'Color',[0 .8  .8])
title(['Rebounding x-velocity';
       '                     '])
xlabel('Time (ms)')
ylabel('Position (mm)')
hold on
plot(t2(1:2:length(t2)),xmaft(1:2:length(t2)),'ow','MarkerSize',4,'MarkerFaceColor','k','LineWidth',0.25)
legend('1^s^t Order Fit','Data','Location','SouthEast')
set(gca, 'box', 'off')
set(gca,'fontsize',20)
hold off

subplot(1,3,3) %vy after detach
plot(t2,y2,'LineWidth',6,'Color',[0 .8  .8])
title(['Rebounding y-velocity';
       '                     '])
xlabel('Time (ms)')
ylabel('Position (mm)')
hold on
plot(t2(1:2:length(t2)),ymaftadj(1:2:length(t2)),'ow','MarkerSize',4,'MarkerFaceColor','k','LineWidth',0.25)
legend('2^n^d Order Fit','Data','Location','SouthEast')
hold off
set(gca, 'box', 'off')
set(gca,'fontsize',20)

prt=1; %print a final plot of the regressions?
if prt==1
    name=['PLOT_RC_VnVt_Model'];
% Defaults for this blog post
width = 20;     % Width in inches
height = 6;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(name,'-dpng','-r300');
end
%---------END Plotting all three measured velocities---------------------

%----------output data for restitution coefficient--------------------
vx2=xfit2(1); %velocity (at t=0) is the second term of the polynomial
vy2=yfit2(1);
vy1=-yfit1(2); %flip to a positive value (drop falling)

if theta>0 & dym~=0 %if surface not flat
%vy1gadj=y-velocity before impact, adjusted for gravity. I include this to
% account for effect of gravity and the change in Y-position before and 
%after impact
vy1gadj=sqrt(vy1^2+2*dym/1000*9.8);
else
vy1gadj=vy1; %no change
end

vn1=vy1gadj*cosd(theta)+vs*sind(theta); %normal velocity before impact
vt1=vs*cosd(theta)-vy1gadj*sind(theta);%tangential velocity before impact
vn2=vx2*sind(theta)+vy2*cosd(theta)-vs*sind(theta);%vn after detach
vt2=vy2*sind(theta)-vx2*cosd(theta)+vs*cosd(theta);%vt after detach
AOI=atand(vt1/vn1); %and calculate angle of incidence

rcn=vn2/vn1; %normal restitution coefficient
rct=vt2/vt1; %tangential RC
rc=sqrt(vn2^2+vt2^2)/sqrt(vn1^2+vt1^2); %overall RC

end %----------------------End regression-------------------------------

disp(['AOI=      ',num2str(AOI)]) %for user
disp(['vn1=      ',num2str(vn1)])
disp(['dY =      ',num2str(dym)])
disp(['RCt=      ',num2str(rct)])
disp(['RCn=      ',num2str(rcn)])
disp(['RC=      ',num2str(rc)])

clipout=[vn1,vy1gadj,AOI,rc,rct,rcn]; %final output vector
disp(kills) %display for user
