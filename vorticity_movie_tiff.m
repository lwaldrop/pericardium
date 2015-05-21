%function F = vorticity_movie
clear all
close all
%These values are smaller than the simulations since not all of the
%vorticity is outputted
length = 0.015; %length of the fluid domain
M = 512;
N = 512;

dx = length/(M);
ds = dx/2;

Ltube = 0.00791; %length of the straight part of the tube
R2 = 0.0024;  %outer radius of the ends of the tube.
Nact = floor(Ltube/ds); %Number of points along the bottom of tube
Nside = floor(pi*R2/ds); %Number of points along the side of the tube

WI = Nact+Nside+Nact+Nside;
L = 360;

%load the vorticity, marker, and forces here
v = load('press_valveless_test');
d = load('markers_valveless_test');
part = load('particles_valveless_test');
p = load('pmarkers_valveless_test');

%x and y are the spatial coordinates
for i = 1:M,
    for j = 1:N,
        x(i,j)=30*(j-1)*(0.0005/(M));
        y(i,j)=30*(i-1)*(0.0005/M);
    end
end

%This determines the maximum value of the vorticity for scaling the plot
cmax = max(v)
cmaxx = max(cmax)

%Create the name of the vorticity file here. Must always make this new or
%it won't output anything.
%aviobj=VideoWriter('womers_movie2.avi');
%open(aviobj);
%'Position',[1 .1 .1 .4],...
%k going from 1 to 20 controls the number of frames in the movie.
%for k=1:150, %150
k =1;
hFig=figure('Visible','Off');
axes1 = axes('Parent',hFig,'YTickLabel',{},'YTick',zeros(1,0),...
    'XTickLabel',{},...
    'XTick',zeros(1,0),...
    'Position',[0.040133779264214 0.0337423312883436 0.905764999095297 0.932515337423313],...
    'FontSize',14,...
    'FontName','Arial');

	%k=1;
	pcolor(x, y, v((k-1)*M+1:(k)*M, 1:N)), shading interp
    colormap(bone)
	caxis([-(cmaxx/16) (cmaxx/16)]); %Use this to control how saturated the vorticity plot is.
	box(axes1,'on');
	hold(axes1,'all');
	hold on
	%plot the wings
	%figure(1)

	plot(d((k-1)*WI+1:(k)*WI,2), d((k-1)*WI+1:(k)*WI,1), '-k', 'LineWidth', 2.0);
	plot(d(((k-1)*WI)+1:k*WI,4), d(((k-1)*WI)+1:k*WI,3), '-k', 'LineWidth', 2.0);
	plot(part(((k-1)*L)+1:k*L,2), part(((k-1)*L)+1:k*L,1), '.w', 'MarkerSize', 6);
	%plot(p(:,2), p(:,1), '-k', 'LineWidth', 1)
	%plot(p(:,4), p(:,3), '-k', 'LineWidth', 1)
	plot(p(1:605,2), p(1:605,1), 'LineWidth', 2,'Color',[0.24705882370472 0.24705882370472 0.24705882370472])
	plot(p(1:605,4), p(1:605,3), 'LineWidth', 2, 'Color',[0.24705882370472 0.24705882370472 0.24705882370472])
	
    hold off
	axis equal;
	axis([0 .015 .002 .01]);
	%xlabel('distance');


	%This part makes the movie.

	%img = hardcopy(hFig, '-dzbuffer', '-r0');
    %writeVideo(aviobj, im2frame(img));
   
    I = getframe(gcf);
    imwrite(I.cdata,['peri_time',num2str(k),'.tif'],'TIFF','Compression','none')
    imwrite(I.cdata,['peri_time',num2str(k),'.png'],'XResolution', 500,'YResolution',500)
    
   close all
   
%end

%close(aviobj)

