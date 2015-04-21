% Analyze data from electro_heart_tube_peri_ciona.c

clear all
%close all


% loads data
a = load('ave_vel_test');
numpts = 27;
mindist=50;

% Calculates number of time steps
n = length(a)/numpts;

% Allocates space
x = zeros(n,numpts);
y = zeros(n,numpts);
peakx = zeros(n,1);
peaky = zeros(n,1);
meany = zeros(n,1);

figure(1)
hold on


be =0;

for i = 1:n
     
    bs = be+1;
    be = bs+(numpts-1);
    x(i,:) = a(bs:be,2);            % Selects the time for this step
    y(i,:) = a(bs:be,1);            % Selects all velocities for step.
    
    [pk,loc] = max(abs(y(i,:)));    % Finds index of peak speed 'loc'
    peaky(i,1) = y(i,loc);          % Uses loc to find speed value.
    peakx(i,1) = max(x(i,:));            % Finds time for current step.
    
    meany(i,1) = mean(y(i,:));
    
    % plots all velocities for this time step on one graph.
    plot(y(i,:))
    
end

hold off

grand_mean = mean(meany)

figure(2)

plot(peakx,peaky)
hold on
plot(peakx,meany,'r-')

[pks,local] = findpeaks(peaky,'minpeakdistance',mindist)
 
plot(peakx(local),peaky(local),'r*')

hold off

tiny_mean=mean(pks)

%State what was found:



% Saves peak data so it can be plotted in R.
fid = fopen('peak_data.csv', 'w') ;
fprintf(fid,'%s\n','peakdata') ;
fprintf(fid, '%s\n') ;
fclose(fid) ;
dlmwrite('peak_data.csv', peaky, '-append') ;

