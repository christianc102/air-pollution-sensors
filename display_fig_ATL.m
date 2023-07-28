function [C, output_args ] = display_fig_ATL( x, mask, pivot, cbounds)
%display_fig Template for displaying image data
%   Detailed explanation goes here
    snapshot = NaN*zeros(113*135,1);
    snapshot(mask==1) = x;
    
    sensors = zeros(113,135);
    P = zeros(size(x)); P(pivot)=ones(size(P(pivot)));
    sensors(mask==1) = P;
    
    C = reshape(real(snapshot),113,135)';
    if (~isempty(cbounds))
        b = imagesc(C,cbounds);
    else 
        b = imagesc(C);%,[-1.799 30.77999]);
    end
    shading interp, colormap(jet), drawnow
    set(b,'AlphaData',(~isnan(C)));
    if (~isempty(sensors))
        hold on
        S = reshape(sensors,113,135)';
        [I,J] = find(S>0);
        
        scatter(J,I,'b.');
        
    end
    axis off
end

