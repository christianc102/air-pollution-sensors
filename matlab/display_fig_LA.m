function [C, output_args ] = display_fig_LA( x, mask, pivot, cbounds)
%display_fig Template for displaying image data
%   Detailed explanation goes here
    snapshot = NaN*zeros(114*86,1);
    snapshot(mask==1) = x;
    
    sensors = zeros(114,86);
    P = zeros(size(x)); P(pivot)=ones(size(P(pivot)));
    sensors(mask==1) = P;
    
    C = reshape(real(snapshot),114,86)';
    if (~isempty(cbounds))
        b = imagesc(C,cbounds);
    else 
        b = imagesc(C);%,[-1.799 30.77999]);
    end
    shading interp, colormap(jet), drawnow
    set(b,'AlphaData',(~isnan(C)));
    if (~isempty(sensors))
        hold on
        S = reshape(sensors,114,86)';
        [I,J] = find(S>0);
        
        scatter(J,I,'b.');
        
    end
    axis off
end

