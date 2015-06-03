% Analyze data from electro_heart_tube_peri_ciona.c

clear all
%close all


% loads data
aveVel = load('ave_vel_test');
pressLine = load('press_line_test');
numpts = 27;
mindist=50;

% Calculates number of time steps
n = length(aveVel)/numpts;

% Allocates space
xVel = zeros(n,numpts);
yVel = zeros(n,numpts);
peakxVel = zeros(n,1);
peakyVel = zeros(n,1);
meanyVel = zeros(n,1);

xPress = zeros(n,numpts);
yPress = zeros(n,numpts);
peakxPress = zeros(n,1);
peakyPress = zeros(n,1);
meanyPress = zeros(n,1);


figure(1)
hold on


be =0;

%Calculations for Average Velocity.
for i = 1:n;
     
    bs = be+1;
    be = bs+(numpts-1);
    xVel(i,:) = aveVel(bs:be,2);            % Selects the time for this step
    yVel(i,:) = aveVel(bs:be,1);            % Selects all velocities for step.
    
    [pk,loc] = max(abs(yVel(i,:)));    % Finds index of peak speed 'loc'
    peakyVel(i,1) = yVel(i,loc);          % Uses loc to find speed value.
    peakxVel(i,1) = max(xVel(i,:));            % Finds time for current step.
    
    meanyVel(i,1) = mean(yVel(i,:));
    
    % plots all velocities for this time step on one graph.
    plot(yVel(i,:))
    
end

hold off

grand_meanVel = mean(meanyVel)

figure(2)

plot(peakxVel,peakyVel)
hold on
plot(peakxVel,meanyVel,'r-')

[pks,local] = findpeaks(peakyVel,'minpeakdistance',mindist);
 
plot(peakxVel(local),peakyVel(local),'r*')

hold off

tiny_meanVel=mean(pks)

%State what was found:



% Saves peak data so it can be plotted in R.
fid = fopen('peak_Veldata.csv', 'w') ;
fprintf(fid,'%s\n','peakVeldata') ;
fprintf(fid, '%s\n') ;
fclose(fid) ;
dlmwrite('peakVel_data.csv', peakyVel, '-append') ;

clear bs be i pks local loc pk

be =0;
figure(3)
hold on

%Calculations for Pressure.
for i = 1:n
     
    bs = be+1;
    be = bs+(numpts-1);
    xPress(i,:) = pressLine(bs:be,2);            % Selects the time for this step
    yPress(i,:) = pressLine(bs:be,1);            % Selects all velocities for step.
    
    [pk,loc] = max(abs(yPress(i,:)));    % Finds index of peak speed 'loc'
    peakyPress(i,1) = yPress(i,loc);          % Uses loc to find speed value.
    peakxPress(i,1) = max(xPress(i,:));            % Finds time for current step.
    
    meanyPress(i,1) = mean(yPress(i,:));
    
    % plots all velocities for this time step on one graph.
    plot(yPress(i,:))
    
end

hold off

grand_meanPress = mean(meanyPress)

figure(4)

plot(peakxPress,peakyPress)
hold on
plot(peakxPress,meanyPress,'r-')

[pksPress,local] = findpeaks(peakyPress,'minpeakdistance',mindist);
 
plot(peakxPress(local),peakyPress(local),'r*')

hold off

tiny_meanPress=mean(pksPress)

%State what was found:



% Saves peak data so it can be plotted in R.
fid = fopen('peakPress_data.csv', 'w') ;
fprintf(fid,'%s\n','peakPressdata') ;
fprintf(fid, '%s\n') ;
fclose(fid) ;
dlmwrite('peakPress_data.csv', peakyPress, '-append') ;

