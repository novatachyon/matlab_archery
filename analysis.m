close all

%% To use: 
%Take a picture of the back of the archery target that you have used. it is
%best if you have the target facing you upright and then you flip it upside
%down towards you. If you flip it sideways then the final figure will be
%rotated (the analysis will be fine, but your target won't be rightside
%up).  Any resolution should be fine, however, the image used by matlab
%needs to be square (you can crop to a 1:1 resolution on any smartphone or
%on your computer). Make sure to try to only have the white paper in the
%cropped picture and not any background or it could edge-detect that
%boundary as shots taken. 

%When the program is run, you should look and see if there are any marks
%near the outer edge of the picture where there were no shots.  If there
%are, you can adjust the tolerance to try to prevent the program from
%counting those; if the problem persists, just crop the image a bit
%smaller. The target size and ring widths are pretty standard but if you
%happen to use a smaller target or something you can adjust these as you
%see fit (Note that the positions of the numbers 10-5 might not line up if
%you change the ring values too much; this can be changed later in the
%program by adjusting the positions of the text pos10, pos9 etc. 

%Note that this program does not know the difference between a large hole
%created by three arrows, and a small hole created by one arrow, because it
%just counts the black spots. To remedy this, you can do two things; a) you
%can place tiny strips of paper over the middle of big holes to split them 
%into two smaller holes that will be counted correctly or b) you can draw 
%little white lines in paint or something to separate bigger holes. There may
%be a way to teach it to do this automatically but I haven't figured it out
%yet!
 
%Note also that usually the image may not be perfectly square i.e. nxn even
%if it has 1:1 dimensions; code in the second section is set up to remove
%the first few rows and a column or two to fix this. If you get an error
%with 'Matrix dimensions must agree' when adding or something, look to the
%right at the variable image and make sure it is nxn. If it is not, you
%can alwasy add in a target(1,:)=[] (which takes out the first row of data)
%or a target(:,1)=[] (which takes out the first column of data) in order to
%make the array square.

%Note that I have assumed a 10/9/8/7/6/5/0 point system with arrows missing
%the target completely not being counted at all.

%The final output of the program is:
%Pointstats: #of shots,  total points, possible points, avg point per shot
%Spreadstats: (all in inches) avg x position, x std dev, avg y pos, y std 
%Quadsats: %shots in I, % in II, % in III, % in IV
%Centering stats: percentage of shots that were 10, 9, 8 ,7, 6, 5, 0 points

%The data is also there to analyze average positions within each ring
%(found in the Ten, Nine, Eight, etc arrays), but it didn't seem useful to
%me.  If you have any questions, my reddit is u/novatachyon!


%% Initial input of the picture and parameters 

file = '2017_4_8.jpg'; %file to be used

tolerance = .2; %sets grayscale tolerance for finding a hole in target:
%should be between 0 and 1, usually lower helps prevent false shots towards
%the edges of the picture being counted.
 
targetsize = 17; %target width in inches
ten = .70; %ten radius and so on, all in inches (all approximate)
nine = 1.6; 
eight = 3.22;
seven = 4.8;
six = 6.35;
five = 8.05;
w = .02; %width of rings to draw on the final figure

image = imread(file);  %reads in image, converts to grayscale and then double
target = rgb2gray(image);
target = double(target);
image =fliplr(rot90(rot90(image))); 

%% Fixing some common issues with the images
for k = 1:8  %takes out the first few rows of data (for some reason, I always have bad edges at the top of my image)
    target(:,1)=[];
    target(length(target),:)=[];
end
target(1,:)=[];
target(:,1)=[];
target(:,1)=[];

target = rot90(target);
target = rot90(target); %rotates the picture if it was flipped upside down
target = fliplr(target); %flips the picture since the picture was taken of the back

target = target./max(max(target)); %normalizes the image and makes the holes white (i.e.=1)
target(target >= tolerance)=2;
target(target < tolerance)=1;
target(target >= 2)=0;

%% Finding the holes and coordinates

holes = edge(target); %finds all areas of holes
c = regionprops(holes,'centroid'); %calculates centroid of each hole
centroids = cat(1,c.Centroid); %concatentates the centroid values

n = length(target); %number of total pixels
conv = targetsize/n; %conversion factor [inches/pixel]
X = linspace(-targetsize/2,targetsize/2,n); %cartesian grid
Y = X;
[xm , ym] = meshgrid(X,Y);
[t , rho] = cart2pol(xm,ym); %converts to a polar grid
m = targetsize/(n-1); %linear spacing of cartesian grid slope
b = -targetsize/2 - m; %makes sure the end/start points are correct
coords = zeros(length(centroids),2); %initializes
for i = 1:length(centroids) %converts centroid pixel values into cartesian values
   coords(i,1) = m*centroids(i,1)+b;
   coords(i,2) = -(m*centroids(i,2)+b);
end

ring1 = zeros(n); ring1(rho<ten+w) =1;ring1(rho<ten-w)=0; %creates arrays for the target rings to be put on the figure
ring2 = zeros(n); ring2(rho<nine+w) =1;ring2(rho<nine-w)=0;
ring3 = zeros(n); ring3(rho<eight+w) =1;ring3(rho<eight-w)=0;
ring4 = zeros(n); ring4(rho<seven+w) =1;ring4(rho<seven-w)=0;
ring5 = zeros(n); ring5(rho<six+w) =1;ring5(rho<six-w)=0;
ring6 = zeros(n); ring6(rho<five+w) =1;ring6(rho<five-w)=0;


%% shot analysis: 

Ten=[];Nine=[];Eight=[];Seven=[];Six=[];Five=[];Zero=[]; %checks to see how many shots were in each point range
%each capitalized array will have the coordinates of a shot that was within
%its ring
for i=1:length(coords)
   x=coords(i,1);
   y=coords(i,2);
   r = sqrt(x^2+y^2);
   if r<=ten
       Ten = [Ten; [x y]];
   elseif r<=nine
       Nine= [Nine; [x y]];
   elseif r<=eight
       Eight= [Eight; [x y]];
   elseif r<=seven
       Seven= [Seven; [x y]];
   elseif r<=six
       Six= [Six; [x y]];
   elseif r<=five
       Five= [Five; [x y]];
   else
       Zero= [Zero; [x y]];
   end
end

cx = mean(centroids(:,1)); %calculates mean/std dev of the centroid values
cy = mean(centroids(:,2));
sigx = std(centroids(:,1));
sigy = std(centroids(:,2));

quad1=0; quad2=0; quad3=0; quad4=0; %calculates how many shots were in each quadrant (I II III IV)
for i =1:length(coords)
    x = coords(i,1); y =coords(i,2);
    if x>0 && y>0
        quad1 = quad1+1;
    elseif x<0 && y>0
        quad2 = quad2+1;
    elseif x<0 && y<0
        quad3=quad3+1;
    else
        quad4 = quad4+1;
    end
end

nshots = length(centroids); %total numbe of shots
totpoints = length(Ten)*10+length(Nine)*9+length(Eight)*8+length(Seven)*7+length(Six)*6+length(Five)*5; %total number of points
posspoints = 10*nshots ; %nshots*10
avgshot = totpoints/nshots; %avg shot point value
avgx = mean(coords(:,1)); %avg coord for all shots and std dev
stdx = std(coords(:,1));
avgy = mean(coords(:,2));
stdy = std(coords(:,2));
percq1 = quad1/nshots*100; %percentage of shots in each quadrant
percq2 = quad2/nshots*100;
percq3 = quad3/nshots*100;
percq4 = quad4/nshots*100;
perctens = length(Ten)/nshots*100; %percentage of shots in each point range
percnines = length(Nine)/nshots*100;
perceights = length(Eight)/nshots*100;
percsevens = length(Seven)/nshots*100;
percsixes = length(Six)/nshots*100;
percfives = length(Five)/nshots*100;
percmiss = length(Zero)/nshots*100;

pointstats = [nshots totpoints posspoints avgshot] %displays the stats as four rows
spreadstats = [avgx stdx avgy stdy]
quadstats = [percq1 percq2 percq3 percq4]
centeringstats = [perctens percnines perceights percsevens percsixes percfives percmiss]


%% Final output plots:

figure(1) %original image of the target back
imshow(image)

figure(2)
fig=target+ring1+ring2+ring3+ring4+ring5+ring6; %combines target image with the rings
str = cell(6,1); %creates a string with the point values and adds the quadrants
for ii=5:10
    str{ii-4} = num2str(ii);
end
str{7}='I';
str{8}='II';
str{9}='III';
str{10}='IV';

f = .4; %correction factor for text
l = length(fig); %calculates relative positions needed for the text overlays (10-5 for rings and I-IV for quadrants);
%note that all of the ring numbers have the same y-position posy

posIa = floor(.94*l); posIb = floor(.04*l); posIIa = posIb; posIIb = posIb; posIIIa = posIb; posIIIb = posIa; posIVa = posIa; posIVb = posIa;
posy = floor(.475*l); pos10 = posy + floor(ten/targetsize*l*f); pos9 = pos10+.07*l; pos8 = pos9 +.09*l; pos7 = pos8+.09*l; pos6 = pos7+.095*l; pos5 = pos6+.1*l;

figfig = insertText(fig,[pos5 posy; pos6 posy; pos7 posy; pos8 posy; pos9 posy; pos10 posy; posIa posIb; posIIa posIIb; posIIIa posIIIb; posIVa posIVb],str,'FontSize',28,'TextColor','white','BoxOpacity',0);
imshow(figfig) % inserts the text strings onto the image and then shows the image

hold on
plot(centroids(:,1),centroids(:,2),'rx','linewidth',2.6)  %plots the centroid positions on the image
plot(cx,cy,'b*','markers',10,'linewidth',2) %plots the average shot on the image
rectangle('Position',[cx-sigx cy-sigy sigx*2 sigy*2],'Curvature',[1 1],'EdgeColor','b','linewidth',1.5) %plots a 1std dev ellipse
rectangle('Position',[cx-2*sigx cy-2*sigy sigx*4 sigy*4],'Curvature',[1 1],'EdgeColor','b','linewidth',1.5) %plots a 2std dev ellipse
hold off


