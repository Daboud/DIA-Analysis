%This program finds the drop's falling velocity in units of pixels/frame,
%and its size, in pixels. This is designed for a STATIONARY surface
%The output will be automatically copied to the user's clipboard, as a vector in with entries:
%[Filename, slant, aspect ratio, Diameter, Velocity];
%The slant is my own arbitrarily scaled measurment of how diagonal the drop
%looks. images with slant above 0.03 should be trashed. The aspect ratio is
%the maximum seen during the video. I only keep AR<0.9
%The output is automatically copied to clipboard, and is formatted to be
%easily pasted into Excel. In Excel, you need cells corresponding to the
%calibration constant of the camer (mm/pixel) and the framerate, to
%generate the dimensional versions of the velocity and diameter

%The program functions using a protocol (bwareopen) that identifies groups of pixels in
%the video, and aims to recognize the drop and not the surface.

clear
exes=[];

f1=1; %First frame to start measuring
% bal=30; %Threshold value, separating black/white in the processed image, for measuring the drop's size and sphericity
blackbar=[]; %This function creates an artificial black bar slashing horizontally through the image. It is helpful when applying a threshold on the image creates an artefact like a floating black spot. Use this to puncture black pixel groups together, so that the program recognizes them as being joined!
sep=1;%Like above, this function helps the program to not misidentify the drop. If there is a black border all the way around the video, set this to =1 to create a white bar thru the top of the image to separate the groups of pixels
balc=60; %threshold value, to find when the drop connects to the surface. Keep about 2x higher then bal
LB = 300; %Lower and upper bounds for pixel group identification. This range should be set to identify the drop and ignore the surface, in terms of their pixel count. Thus in the videos, the surface should be much large than the area of the drop so the separation is easily recognized.
UB = 2000; 
fskip=50; %how often to refresh plot to view progress
Filename=['H200 60o w9 07'];

use=[Filename '.avi'];
v = VideoReader(use);
frames=round(v.Duration*v.FrameRate);
f2=round(frames/3); %this attempts to estimate how much of the video to load (to save time). Change the 3 to other values if your video is cut short, or if it takes a long time to load
vid = read(v,[f1 f2]);


for i=[f1:f2]%main loop
    disp(i)
    I0=vid(:,:,:,i-f1+1);
    [m,n]=size(I0); %height, width of image
    if i==1; %This set of conditions sets a threshold brightness value based on the ovrall brightness of the video. You may need to change these if your videos are much darker or brighter than mine. The goal is for the drop diameter measurement to not be affected by the brightness.
        brightn=max(max(I0));
        if 235>brightn & brightn>225
            bal=25
        elseif brightn<=225
            bal=30
        elseif brightn>235
            bal=20
        end
    end

    I=I0; %Make a black/white version of the image using the threshold value
    black=find(I<bal); I(black)=0;
    white=find(I>=bal); I(white)=255; %turn black & white

    Ic=I0; %this black/white version is used to identify when the drop contacts the surface
    black=find(Ic<balc); Ic(black)=0;
    white=find(Ic>=balc); Ic(white)=255; %turn black & white

    if sep==1 %creates a white separation line on top.
        Ic(1:200,round(n/1.7))=255;
        I(1:200,round(n/1.7))=255;
    end

    I=255-I; I=imfill(I); I=255-I; %Image Filling, to remove white spot in drop
    Ic=255-Ic; Ic=imfill(Ic); Ic=255-Ic; %Image Filling

    I(1:m,end)=0; I(end,1:n)=0; I(1:m,1)=0; I(1,1:n)=0; %This adds a black border, to try and bridge together all elements that are not the drop we are searching for
    I(blackbar,:)=0;
    Ic(1:m,end)=0; Ic(end,1:n)=0; Ic(1:m,1)=0; Ic(1,1:n)=0; 
    Ic([blackbar,blackbar+1],:)=0;

    Ibin=(255-I)/255; %convert I to binary information
    Iout = xor(bwareaopen(Ibin,LB),  bwareaopen(Ibin,UB)); %bwareopen is a useful function that identifies all groups of pixels, as well as their pixel count and indexed locations. This is the most important function in the code
    Id=(1-Iout)*255; %create image using only kept pixel group (hopefully the drop)

    [y,x]=find(Id==0); %find black pixels

    Ibinc=(255-Ic)/255; %convert I to binary information
    Ioutc = xor(bwareaopen(Ibinc,LB),  bwareaopen(Ibinc,UB)); %keep only between
    Idc=(1-Ioutc)*255; %create image using only kept (drops)

    [yc,xc]=find(Idc==0); %find black pixels

    if numel(xc)==0 %this condition identifies when the drop hits the surface, ending the main loop
        f2=i-1;
        break
    end

    xmp(i-f1+1)=mean(x); %get mean coordinates of drop location
    ymp(i-f1+1)=mean(y);
    exes=[exes;xmp(i-f1+1) ymp(i-f1+1)]; %matrix is used for center of mass

    %-----------------  Start AR Measurement  -----------------------
    %This section measures the aspect ratio of the drop
    [yy,xx]=find(Id==0); %find all black pixels
    xleft=min(xx); %find leftmost position...
    xright=max(xx);
    yroof=min(yy);
    yfloor=max(yy);

    dx=xright-xleft+1; %width
    dy=yfloor-yroof+1; %height

    Rgirth=length(find(xx==xright)); %How many pixels constitute the edge?
    Lgirth=length(find(xx==xleft));
    Tgirth=length(find(yy==yroof));
    Bgirth=length(find(yy==yfloor));

    if Rgirth<3 | Lgirth<3 %If 2 or fewer, truncate by half a pixel
        dx=dx-0.5;
    elseif Bgirth<3 | Tgirth<3
        dy=dy-0.5;
    end

    if dy/dx<=1 %measure aspect ratio of drop
        AR(i-f1+1)=dy/dx;
    else
        AR(i-f1+1)=dx/dy;
    end
    %-----------------  End AR Measurement  -----------------------

    %-----------------  Start Slant Measurement  -----------------------
    numlines=floor(round(dy)/2); %how many rows tall is the top and the bottom of the drop?

    roofhalf=find(yy<=yroof+numlines-1); %find where in the yy vector the roof pixels are
    roofpixy=yy(roofhalf); %identify in x-y coordinates the y positions
    roofpixx=xx(roofhalf);% and the x positions
    roofmidx=mean(roofpixx); %average them to find midpoint of roof
    roofmidy=mean(roofpixy);

    floorhalf=find(yy>=yfloor-numlines+1); %find where in the yy vector the floor pixels are
    floorpixy=yy(floorhalf); %identify in x-y coordinates the y positions
    floorpixx=xx(floorhalf);% and the x positions
    floormidx=mean(floorpixx); %average them to find midpoint of floor
    floormidy=mean(floorpixy);

    obl(i-f1+1)=abs(floormidx-roofmidx); %pixel difference between L and R
    %-----------------  End Slant Measurement  -----------------------

    D(i-f1+1)=sqrt(dx*dy); %measured diameter of drop

    if i==f1   |   i==f2  | ((i)/fskip==round(i/fskip)) %pause regularly to show progress

        colormap(gray);
        MainFig = figure;
        set(MainFig, 'Position', [200, 50, 1200, 350]);
        close(1) 

        subplot(1,3,1)
        imshow(I0) %original image
        hold on
        subplot(1,3,2)
        image(Ic) %threshold applied image, seeking drop-surface contact
        hold on
        subplot(1,3,3)
        image(Id) %should be only drop
        colormap(gray)
        hold on
        plot(roofmidx,roofmidy,'or','MarkerSize',5,'LineWidth',1)
        plot(floormidx,floormidy,'or','MarkerSize',5,'LineWidth',1)

        subplot(1,3,3)
        hold on
        plot(xmp(end),ymp(end),'xb','MarkerSize',15,'LineWidth',3)    

        pause

    end
end

image(Ic) %when loop breaks, show final frame to confirm contact properly identified
pause
close

Tpix=[f1:f2]; %frame timeline
D(find(D==0))=[];

Pd=polyfit(Tpix,ymp,2); %second order polynomial fit to drop's y-location during fall
modd=polyval(Pd,Tpix); %modeled position line
Vd=Pd(1)*2*Tpix(end)+Pd(2); %slope of modeled line at moment of contact
dx1=1; dx2=Tpix(end); %tangent line from final point showing slope
dy2=ymp(end); dy1=dy2-Vd*(dx2-dx1);

Du=unique(D); %find each unique Diameter found
lu=length(Du);
if lu==3 %This procedure selects a measurement of D that best describes the overall set of measurements
    choice=2
end
if lu==2
    choice=1
end
if lu==1
    choice=1
end
if lu>3
    choice=round(0.15*length(Du))+1; %pick one ~20% of the way up
end
Dgood=Du(choice); 
goods=find(D==Dgood);
reds=ones(1,length(goods))*Dgood; %highlight in red that value

MainFig = figure; %final plotting
set(MainFig, 'Position', [200, 50, 1300, 400])
subplot(1,3,1)
plot(Tpix,modd,'-c','linewidth',5)
hold on
plot(Tpix,ymp,'.k','markersize',7)
plot([dx1 dx2],[dy1 dy2],'-k')
title('Drop Velocity')
subplot(1,3,2)
plot(Tpix,D,'ok')
hold on
plot(goods,reds,'or','markersize',10,'linewidth',2)
plot([1 f2], [Dgood Dgood],'-r')
title('Drop Diameter (pix)')
subplot(1,3,3)
yyaxis left
plot(Tpix,AR,'-b')
ylabel('Aspect Ratio')
ylim([0 1])
yyaxis right
plot(Tpix,obl/Dgood,'-r')
hold on
ylabel('Slant')
ylim([0 0.1])

clipout=[Filename, '	', num2str(max(obl/Dgood)),'	', num2str(min(AR)),'	',num2str(Dgood),'	', num2str(Vd)];
clipboard('copy',num2str(clipout)) %copy to clipboard vector with measurements
