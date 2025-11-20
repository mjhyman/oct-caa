%% scale bar function
function sbar(sbar_len, sbar_thick, vox, im)
% Add scale bar to bottom right corner
% INPUTS:
%   - sbar_len (uint): scale bar length (microns)
%   - sbar_thick (uint): scale bar line width (builtin linewidth)
%   - vox (vector): isotropic voxel size (microns)
%   - slice (float matrix): single depth of image to display

%%% Define size
% scalebar length in pixels
sbar_px = sbar_len ./ vox;
% get image size
[imHeight, imWidth] = size(im,[1,2]);

%%%% Position: bottom right margin
% small margin (2% of width)
x_end = imWidth - round(imWidth*0.02);
x_start = x_end - sbar_px;          
% a little above bottom (3% of height)
y_pos = imHeight - round(imHeight*0.03);

%%% Draw scale bar (white line)
hold on;
plot([x_start x_end], [y_pos y_pos], 'w', 'LineWidth', sbar_thick);
hold off;
% Disable x,y ticks
set(gca,'XTick',[]); set(gca,'YTick',[])

end