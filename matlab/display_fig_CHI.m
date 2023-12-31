function [C, output_args ] = display_fig_CHI( x, mask, pivot, cbounds)
%display_fig Template for displaying image data
%   Detailed explanation goes here
    snapshot = NaN*zeros(137*130,1);
    snapshot(mask==1) = x;
    
    sensors = zeros(137,130);
    P = zeros(size(x)); P(pivot)=ones(size(P(pivot)));
    sensors(mask==1) = P;
    
    C = reshape(real(snapshot),137,130)';
    if (~isempty(cbounds))
        b = imagesc(C,cbounds);
    else 
        b = imagesc(C);%,[-1.799 30.77999]);
    end
    shading interp, colormap(jet), drawnow
    set(b,'AlphaData',(~isnan(C)));
    if (~isempty(sensors))
        hold on
        S = reshape(sensors,137,130)';
        [I,J] = find(S>0);
        
        scatter(J,I,'b.');
        
    end
    axis off
end

