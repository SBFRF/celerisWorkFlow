function colorbar_image = create_colorbar_jpg(cAxis)
%KML.COLORBAR(cAxis,colorMap) Create a color bar visible as a screen overlay.
%   Similar to built-in colormap function, input arguments shoud be the 
%   minimum and maximum range values, and colorMap a matrix such as the one
%   created by calling the built-in function JET
%
%   Copyright 2013 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.6 $  $Date: 2013/05/17 17:17:17 $
ss = get(0,'ScreenSize')*.3;


l = linspace(cAxis(1),cAxis(2),1000);
fh = figure;


ratio = 100;

pcolor([0*l.' 1+0*l.'],[l.' l.'],[l.' l.'])
ylabel('Ocean Elevation (m)')

set(gcf,'Position',[0 0 ss(4)./ratio  1.5*ss(4)]);%[ss(1)+ss(4)./ratio  ss(2) ss(4)./ratio  ss(4)]);
set(gca,'XTick',[]);
sz         = [2/ratio 1];
overlayPos = [2 1];
screenPos  = [1 1];

movegui(fh,'north');

shading flat;
caxis(cAxis);
colormap(jet);


name = 'Colorbar';
name = [name '.png'];

bgColor = [0 0 0];
fgColor = [1 1 1];

set(gca,'Color',bgColor,'XColor',fgColor,'YColor',fgColor,'YAxisLocation','right');
ytickangle(gca,90);
set(fh,'Color',bgColor);
set(gca,'FontWeight','bold','FontSize',10);

set(fh,'PaperPositionMode','auto','InvertHardcopy','off')
print('-dpng','-r0',name,fh)
colorbar_image = imread(name);

close(fh)