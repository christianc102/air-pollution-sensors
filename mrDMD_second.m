%% Script to display particular modes of an mrDMD map produced from mrDMD_plot script

clear; close all; clc
ncdfpath = 'NetCDFs/';
maskpath = 'Masks/';
figpath = 'Figures/';

% config
urb_area = 'LosAngelesLongBeachAnaheimCA';
pollutant = 'PM';
levels = 4;
max_cycles = 1;

% Set PRINT_FIG=true to export figures
PRINT_FIG = true;
set(0, 'defaultfigurecolor', 'w');

% Read in the data
dat = ncread([ncdfpath, urb_area, pollutant, '.nc'], [pollutant, ' concentration']); 
mask = ncread([maskpath, urb_area, 'mask.nc'] , 'Urban Area');
mask(isnan(mask))=0;
dat(isnan(dat))=0;

% Second 4096 days of data, excluding corrupted dates
dat = dat(:,:,setdiff((2082:6180),[3291,5689,5690]));

%numyears = 16;
%nweeks = numyears*52;


Y = zeros(length(mask(mask==1)),size(dat,3));
for i=1:size(dat,3)
    Band = dat(:,:,i);
    Y(:,i) = Band(mask==1);
end

[m,n] = size(mask);
N = m*n;

% Construct cost function as implemented in sensors-cost-paper code
%f2 = ones(m,n);
%for j = 3:357
%    for jj = 3:177
%        if mask(j,jj) == 0
%            f2(j-2:j+2,jj-2:jj+2) = 0;
%        end
%    end
%end

%race cost fn:
%f2 = ncread([datpath,'PropNWBostonx.nc'],'Band1'); %might need to make this inverse
%f2 = ncread([datpath,'NormalizedIncomeBostonx.nc'],'Band1'); %might need to make this inverse
%f2(isnan(f2))=0;
%f = reshape(f2,N,1);

% mrDMD code has no land pixels, so we remove those from the cost function
%f = f(mask == 1);


[N,M] = size(Y);

% mrDMD frequencies should be yearly
% scale omegas by 52
% but keep dt=1 week.
dt = 1;%mean(diff(time))/7; 
%timeval = time(nweeks+1:end);

% 16 year period ending exactly in 2017
%time2 = time(end-8-831:end-8);
%tree2 = mrDMD_fb(Y(:,end-8-831:end-8),dt,10,1,4, true);

% mrDMD, 16 year period beginning in 1990
%time = setdiff((1:100),[90])';
time = setdiff((2082:6180),[3291,5689,5690]);

%tree = mrDMD_fb(Y(:,1:length(time)),dt,8,1,10, true); %change mrDMD modal tree here
tree = mrDMD_fb(Y(:,1:length(time)), dt, 8, max_cycles, levels, true); %change mrDMD modal tree here

[U,S,V] = svd(Y(:,1:length(time)),'econ');
[~,~,piv] = qr(U(:,1:30)',0);
qdeim = piv(1:30);
%% plot amplitudes of 1990 mrDMD

%indt below adds lines to the tree figure to see when in time series it
%lies. Comment out for final fig. 
%mkelp: create test set of 10%, seed it. See old error code?

%indt = [12, 120, 1000, 300, 900, 100, 2000, 3000, 2500, 3500]; % measure in summer months; mkelp adjust here for test/train set
rng(10)
%indt = randperm(2082, 400) %extrapolate
indt = randperm(4000, 400); %interpolate

[ptree, map, low_f_cutoff] = mrDMD_map(tree);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree{L,J}.T;

T = datetime(2000,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*365);
%set(gca,'XTick', (0:J-1) + 0.5,'XTickLabel',num2str([month(T),year(T)],'%d/%d'));
axis xy; axis tight; box on; grid on

ylim = get(gca,'ylim');

times = year(datetime(2000,1,1,0,0,0) + days(time(indt)));
hold on
%for i=1:length(indt)
%    stem(1+indt(i)/shortInt, ylim(2),'k.','LineWidth',1.5)
%    
%end

set(gca, 'XTickLabelRotation',45,'ylim',ylim,'FontSize',12);
colormap(gca,'pink');shading flat

if PRINT_FIG
    %file_name = strcat(figpath, 'FIG_MRDMD_MAP.fig');
    %savefig(file_name);
    print(strcat(figpath, 'FIG_MRDMD_MAP_', urb_area, pollutant, 'L=', string(levels), 'cyc=', string(max_cycles)), '-dpdf', '-fillpage');
end
%% collect and display unique mrDMD modes

lmap = []; jmap = [];
Lib = []; Lib2 = [];
Omega = []; Omega2 = [];
Amp = []; Amp2 = [];
Periods = [];
Lambda = []; Lambda2 = [];

tol = 1e-2; % 1e-2

%% Plot background mode (l=1,j=1) and save image
for l=1
    for j=1
        Lib = [Lib tree{l,j}.Phi(1:N,:)];        
        Omega = [Omega; tree{l,j}.omega*365]; % yearly (!)
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
        Lambda = [Lambda; tree{l,j}.lambda];

        figure;
        display_fig_LA(abs(tree{l,j}.Phi(1:N,1)),mask,[],[]);
        colorbar;
        if PRINT_FIG
            %file_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.fig');
            %savefig(file_name);
            png_name = strcat(figpath, 'FIG_MRDMD_', urb_area, pollutant, 'BG_MODE_L=',string(l),'_J=',string(j), 'L=', string(levels), 'cyc=', string(max_cycles), '.png');
            saveas(gcf,png_name);
        end
    end
end




%% collect and display unique mrDMD modes

lmap = []; jmap = [];
Lib = []; Lib2 = [];
Omega = []; Omega2 = [];
Amp = []; Amp2 = [];
Periods = [];
Lambda = []; Lambda2 = [];

tol = 1e-2; % 1e-2

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];        
        Omega = [Omega; tree{l,j}.omega*365]; % yearly (!)
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
        Lambda = [Lambda; tree{l,j}.lambda];
        
        if imag(tree{l,j}.omega*365)>tol
            disp('big mode!')
            figure;
            display_fig_LA(abs(tree{l,j}.Phi(1:N,ind(kk))),mask,[],[]);
            colorbar; 
            clim([0, 0.02]);
            str = ['Phi97' num2str(l) ',' num2str(j) '_' ...
                num2str(round(imag(tree{l,j}.omega(ind(kk))*365),2))];
            if PRINT_FIG
                %file_name = strcat(figpath, 'FIG_MRDMD_MODE_L=',string(l),'_J=',string(j),'.fig');
                %savefig(file_name);
                png_name = strcat(figpath, 'FIG_MRDMD_', urb_area, pollutant, '_MODE_L=',string(l),'_J=',string(j), 'L=', string(levels), 'cyc=', string(max_cycles), '.png');
                saveas(gcf,png_name);
            end
        end
    end
end
%       Code from mrDMD paper to produce maps of the modes; however, the 
%       tol value may not be calibrated correctly, so a lot of modes are
%       printed and the program takes too long to run. Use
%       mrDMD_specific_modes_second_half.m
%