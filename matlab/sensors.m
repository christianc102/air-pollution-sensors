%% Script to display particular modes of an mrDMD map produced from mrDMD_plot script

clear; close all; clc
ncdfpath = 'NetCDFs/';
maskpath = 'Masks/';
figpath = 'Figures/Sec/';

% config
urb_area = 'LosAngelesLongBeachAnaheimCA';
pollutant = 'O3NO2ratio';
max_cycles = 1;
levels = 7;

% Set PRINT_FIG=true to export figures
PRINT_FIG = true;
set(0, 'defaultfigurecolor', 'w');

% Read in the data
dat = ncread([ncdfpath, urb_area, pollutant, '.nc'], [pollutant, ' concentration']); 
mask = ncread([maskpath, urb_area, 'mask.nc'] , 'Urban Area');
mask(isnan(mask))=0;
dat(isnan(dat))=0;

% Short time span test
% dat = dat(:,:,(129:256));

% First 4096 days of data, excluding corrupted dates
%dat = dat(:,:,setdiff((1:4097),3291));

% Second 4096 days of data, excluding corrupted dates
dat = dat(:,:,setdiff((2112:6210),[3291,5689,5690]));

%numyears = 16;
%nweeks = numyears*52;


Y = zeros(length(mask(mask==1)),size(dat,3));
for i=1:size(dat,3)
    Band = dat(:,:,i);
    Y(:,i) = Band(mask==1);
end

%[m,n] = size(mask);

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
%time = (129:256);
%time = setdiff((1:4097), 3291);
time = setdiff((2112:6210),[3291,5689,5690]);

[U,S,V] = svd(Y(:,1:length(time)),'econ');
[~,~,piv] = qr(U(:,1:30)',0);
qdeim = piv(1:30);

% finding optimal truncation (r) using optimal SVHT
sigs = diag(S);
beta = size(Y(:,1:length(time)), 1) / size(Y(:,1:length(time)),2);
thresh = optimal_SVHT_coef(beta,0) * median(sigs);
disp(length(sigs(sigs>thresh)))
figure
semilogy(sigs, '-ok', 'LineWidth', 1.5)
grid on
hold on
xlim([0 length(sigs)])
ylim([1 10^6])

semilogy(sigs(sigs>thresh), 'bo', 'LineWidth', 1.5)
plot([-20 length(sigs)], [thresh thresh], 'b--', 'LineWidth', 2)

%tree = mrDMD_fb(Y(:,1:length(time)),dt,8,1,10, true); %change mrDMD modal tree here
tree = mrDMD_fb(Y(:,1:length(time)), dt, length(sigs(sigs>thresh)), max_cycles, levels, true); %change mrDMD modal tree here

% if PRINT_FIG
%     png_name = strcat(figpath, 'optimalSVHT', urb_area, pollutant, '.png');
%     saveas(gcf,png_name);
% end
% results - r 
% first
% LA NO2 = 976
% LA O3 = 877
% LA PM = 1064
% CHI NO2 = 

% second
% LA NO2 = 992
% LA O3 = 911
% LA PM = 1023
% ATL NO2 = 912
% ATL O3 = 992
% ATL PM = 937
% plot amplitudes of 1990 mrDMD

%indt below adds lines to the tree figure to see when in time series it
%lies. Comment out for final fig. 
%mkelp: create test set of 10%, seed it. See old error code?

%indt = [12, 120, 1000, 300, 900, 100, 2000, 3000, 2500, 3500]; % measure in summer months; mkelp adjust here for test/train set
rng(10)
%indt = randperm(2082, 400) %extrapolate
indt = randperm(4000, 400); %interpolate
figure
[ptree, map, low_f_cutoff] = mrDMD_map(tree);
[L, J] = size(tree);

imagesc(-sqrt(map));

shortInt = tree{L,J}.T;

T = datetime(2000,1,1,0,0,0) + days(time(1:shortInt:end));
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', low_f_cutoff*365);
%set(gca,'XTick', (0:J-1) + 0.5,'XTickLabel',num2str([month(T),year(T)],'%d/%d'));
axis xy; axis tight; box on; grid on

ylim = get(gca,'ylim');
% 
% times = year(datetime(2000,1,1,0,0,0) + days(time(indt)));
% hold on
%for i=1:length(indt)
%    stem(1+indt(i)/shortInt, ylim(2),'k.','LineWidth',1.5)
%    
%end

set(gca, 'XTickLabelRotation',45,'ylim',ylim,'FontSize',12);
colormap(gca,'pink');shading flat
colorbar;

if PRINT_FIG
    %file_name = strcat(figpath, 'FIG_MRDMD_MAP.fig');
    %savefig(file_name);
    %print(strcat(figpath, 'FIG_MRDMD_MAP_', urb_area, pollutant, 'cyc=', string(max_cycles)), '-dpdf', '-fillpage');
    png_name = strcat(figpath, 'FIG_MRDMD_SENS_SEC_CALENDAR', urb_area, pollutant, 'cyc=', string(max_cycles), 'lvl=', string(levels), '.png');
    saveas(gcf,png_name);
end

%% collect unique mrDMD modes

lmap = []; jmap = [];
Lib = [];
Omega = [];
Amp = [];
Periods = [];

for l=1:L
    for j=1:2^(l-1)
        Lib = [Lib tree{l,j}.Phi(1:N,:)];
        Omega = [Omega; tree{l,j}.omega*12]; % weekly
        Amp = [Amp; tree{l,j}.P];
        Periods = [Periods; tree{l,j}.T];
    end
end

% Different values of gamma, the cost-balance coefficient. Play with this 
Gammas = [0];  %need to play with this, very sensitive 

%Gammas = [0.000001 0.00001 0.0001 0.001 0.01 0.1 1 2];  %need to play with this, very sensitive 
%Gammas = [0 0.0001 0.001 0.005 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2];  %need to play with this, very sensitive 


%Gammas = [2];  %need to play with this, very sensitive 
%Gammas = [0:10:100];
%Gammas = [0 0.002 0.02 0.04]; %play around with ranges to see what is appropriate
%look at map print out to see a gradual shifting to less white areas

Costs = []; % Stores total cost of sensors for each Gamma value
Errors = []; 

% Invert f; store the original f in new var fcopy
%fcopy = f;
%f(:) = 0.9729./f; %~f; %maybe remove this

for Gamma=Gammas
    [~,~,piv] = qr(Lib.','vector'); % calling qr function in original mrDMD code
    %[~,~,piv] = qrpc(Lib.',Gamma*f); % calling qrpc function from sensors-cost-paper
    sens = piv(1:size(Lib,2)); % store the locations of top size(Lib,2) sensors == mkelp: these are ideal # sensors
    %sens = piv(1:1000);
    %Cost = sum(fcopy(sens)); % Sum of cost function’s value at each of chosen sensors
    %Costs = [Costs Cost]; % Append cost for each Gamma value to Costs
    
    fprintf('ideal number of sensors: %d', size(Lib, 2))
    
    Phi = cell(1,L);
    omega = cell(L,1);
    b = cell(L,1);
    t0 = cell(L,1);

    PhiJ = cell(1,J); win =1;
    lab = [];
    xL = zeros(N,L); %prediction at each level
    update = false;
    Lev_coeff = [];

    Xrecon = [];

    %% compute reconstruction of SST from multiscale sensor measurements
    for i=1:length(time)
        x = Y(:,i);

        %xavg = L(:,i); % filtered data for sensors

        if (i>2 && i < length(time)-1)
            xavg = mean(Y(:,i-2:i+2),2);
        else 
            xavg = x;
        end


        % at each level update t-t0
        for l=1:L
            Int = tree{l,1}.T;

            if (mod(i,Int)==0 && i< length(time))
                t0{l} = time(i+1);
                ind = floor((i+1)/Int)+1;
                update = true;
            elseif (i==1)
                ind = 1;
                t0{l} = time(1);
                update = true;
            end

            if (update)   
                %disp([i ind]) % update modes when time windows change

                Phi{l} = tree{l,ind}.Phi(1:N,:);
                omega{l} = tree{l,ind}.omega;
                b{l} = tree{l,ind}.P;
                %sub{l} = ceil(1/8/pi/tree{l,ind}.rho/dt);


                Modes = cell2mat(Phi);

                ti = (time(i)-t0{l}); 
                xL(:,l) = Phi{l}*(exp(2*pi*omega{l}*ti).*b{l});
                if (i==1)
                    x0 = x;
                else
                    x0 = Y(:,i+1);
                end
                lev_coeff = xL\x0;           

                if (l==L) % store J sets of modes by last window
                    PhiJ{win} = cell2mat(Phi);
                    lab = [lab; win*ones(size(PhiJ{win},2),1)];
                    Lev_coeff = [Lev_coeff lev_coeff];
                    win = win+1;
                end

                update = false;
                r = size(Modes,2);              
                disp([size(Modes,2) length(sens)]);
            end                

            ti = (time(i)-t0{l}); 
            xL(:,l) = Phi{l}*(exp(2*pi*omega{l}*ti).*b{l});
        end


        % weighted mrDMD reconstruction
        % least-squares fit of by-level approximation to true snapshot
        xpred = xL*lev_coeff;    

        % collect Phi at all levels and approximate
        xhat = Modes*(Modes(sens,:)\x(sens));
        xhat2 = Lib*(Lib(sens,:)\x(sens));  

        xpod = U(:,1:r)*(U(qdeim,1:r)\x(qdeim));

        if (ismember(i, indt))
            Xrecon = [Xrecon xhat];
        end
        errpred(i) = norm(x-xpred)/norm(x);
        errsens(i) = norm(x-xhat)/norm(x);
        errsens2(i) = norm(x-xhat2)/norm(x);
        errpod(i) = norm(x-xpod)/norm(x);

    end
    
    % Calculating RMSE and MPE errors using functions from summer PM2.5 code
    rmse_errors = []; 
    mpe_errors = [];
    for i=1:size(Xrecon,2)
        [rmse,mpe] = rmse_mpe(Xrecon(:,i),dat(:,:,indt(i)),mask);
        rmse_errors = [rmse_errors; rmse];
        mpe_errors = [mpe_errors; mpe];
    end
    
    % rmse_errors contain rmse errors at all locations, we take the mean
    rmse_mean = mean(rmse_errors);
    
    % error calculation from original mrDMD code
    errpred_mean = mean(errpred);
    errsens_mean = mean(errsens);
    errsens2_mean = mean(errsens2);
    errpod_mean = mean(errpod);
    
    % Append rmse error for each gamma value to Errors
    Errors = [Errors rmse_mean]; %add different error here if you want
    
    %mkelp: prints a sensor map figure for each gamma 
    figure;
    display_sensors_LA(mask,sens(101:200),'r');
    if PRINT_FIG
        %file_name = strcat(figpath, 'FIG_MRDMD_MAP.fig');
        %savefig(file_name);
        %print(strcat(figpath, 'FIG_MRDMD_MAP_', urb_area, pollutant, 'cyc=', string(max_cycles)), '-dpdf', '-fillpage');
        png_name = strcat(figpath, 'FIG_MRDMD_SENS_SEC_101-200', urb_area, pollutant, 'cyc=', string(max_cycles), 'lvl=', string(levels), 'sens=', num2str(size(Lib, 2)), '.png');
        saveas(gcf,png_name);
    end


    %figure;
    %display_sensors(mask,sens(1:1000),'r');
    %if PRINT_FIG
    %    filename = ['figures/final2/recon_sensors_top1000_' num2str(Gamma)];
    %    print(filename, '-dpdf', '-fillpage');
    %end

    %figure;
    %display_sensors(mask,sens(1:500),'r');
    %if PRINT_FIG
    %    filename = ['figures/final2/recon_sensors_top500_' num2str(Gamma)];
    %    print(filename, '-dpdf', '-fillpage');
    %end

    %figure;
    %display_sensors(mask,sens(1:750),'r');
    %if PRINT_FIG
    %    filename = ['figures/final2/recon_sensors_top750_' num2str(Gamma)];
    %    print(filename, '-dpdf', '-fillpage');
    %end

%     figure;
%     display_sensors_MAD(mask,sens,'r');
%     if PRINT_FIG
%         filename = ['figures/final_20220920_race_STL/recon_sensors_top250_' num2str(Gamma)];
%         print(filename, '-dpdf', '-fillpage');
%     end


%     figure;
%     display_sensors_MAD(mask,sens(1:200),'r');
%     if PRINT_FIG
%         filename = ['figures/final_20220920_race_STL/recon_sensors_top200_' num2str(Gamma)];
%         print(filename, '-dpdf', '-fillpage');
%     end
% 
% 
%     figure;
%     display_sensors_MAD(mask,sens(1:150),'r');
%     if PRINT_FIG
%         filename = ['figures/final_20220920_race_STL/recon_sensors_top150_' num2str(Gamma)];
%         print(filename, '-dpdf', '-fillpage');
%     end
% 
%     figure;
%     display_sensors_MAD(mask,sens(1:100),'r');
%     if PRINT_FIG
%         filename = ['figures/final_20220920_race_STL/recon_sensors_top100_' num2str(Gamma)];
%         print(filename, '-dpdf', '-fillpage');
%     end
% 
% %    sensor_matrix = display_sensors_matrix(mask,sens,'r'); 
% %    %imagesc(sensor_matrix)
% % 
% %    numrow = 101;
% %    numcol = 143;
% %    filename = ['figures/final2/netcdf/netcdf_ideal_' num2str(Gamma) '.nc'];
% %    ncid = netcdf.create(filename,'NC_WRITE');
% %    dimidrow = netcdf.defDim(ncid,'rows',numrow);
% %    dimidcol = netcdf.defDim(ncid,'length',numcol);
% %    varid = netcdf.defVar(ncid,'loc_Boston','NC_DOUBLE',[dimidrow dimidcol]);
% %    netcdf.endDef(ncid);
% %    netcdf.putVar(ncid,varid,sensor_matrix);
% %    netcdf.close(ncid);
% %    ncid2 = netcdf.open(filename,'NC_NOWRITE');
% %    data_copy = netcdf.getVar(ncid2,0);
% %    if isequal(sensor_matrix,data_copy)
% %        disp('Data match');
% %    else
% %        disp('Data mis-match');
% %    end
% %
% %
% %    sensor_matrix = display_sensors_matrix(mask,sens(1:1000),'r');
% %    %imagesc(sensor_matrix)
% %
% %    numrow = 101;
% %    numcol = 143;
% %    filename = ['figures/final2/netcdf/netcdf_top1000_' num2str(Gamma) '.nc'];
% %    ncid = netcdf.create(filename,'NC_WRITE');
% %    dimidrow = netcdf.defDim(ncid,'rows',numrow);
% %    dimidcol = netcdf.defDim(ncid,'length',numcol);
% %    varid = netcdf.defVar(ncid,'loc_Boston','NC_DOUBLE',[dimidrow dimidcol]);
% %    netcdf.endDef(ncid);
% %    netcdf.putVar(ncid,varid,sensor_matrix);
% %    netcdf.close(ncid);
% %    ncid2 = netcdf.open(filename,'NC_NOWRITE');
% %    data_copy = netcdf.getVar(ncid2,0);
% %    if isequal(sensor_matrix,data_copy)
% %        disp('Data match');
% %    else
% %        disp('Data mis-match');
% %    end
% %
% %    sensor_matrix = display_sensors_matrix(mask,sens(1:750),'r');
% %    %imagesc(sensor_matrix)
% %
% %    numrow = 101;
% %    numcol = 143;
% %    filename = ['figures/final2/netcdf/netcdf_top750_' num2str(Gamma) '.nc'];
% %    ncid = netcdf.create(filename,'NC_WRITE');
% %    dimidrow = netcdf.defDim(ncid,'rows',numrow);
% %    dimidcol = netcdf.defDim(ncid,'length',numcol);
% %    varid = netcdf.defVar(ncid,'loc_Boston','NC_DOUBLE',[dimidrow dimidcol]);
% %    netcdf.endDef(ncid);
% %    netcdf.putVar(ncid,varid,sensor_matrix);
% %    netcdf.close(ncid);
% %    ncid2 = netcdf.open(filename,'NC_NOWRITE');
% %    data_copy = netcdf.getVar(ncid2,0);
% %    if isequal(sensor_matrix,data_copy)
% %        disp('Data match');
% %    else
% %        disp('Data mis-match');
% %    end
% %
% %    sensor_matrix = display_sensors_matrix(mask,sens(1:500),'r');
% %    %imagesc(sensor_matrix)
% %
% %    numrow = 101;
% %    numcol = 143;
% %    filename = ['figures/final2/netcdf/netcdf_top500_' num2str(Gamma) '.nc'];
% %    ncid = netcdf.create(filename,'NC_WRITE');
% %    dimidrow = netcdf.defDim(ncid,'rows',numrow);
% %    dimidcol = netcdf.defDim(ncid,'length',numcol);
% %    varid = netcdf.defVar(ncid,'loc_Boston','NC_DOUBLE',[dimidrow dimidcol]);
% %    netcdf.endDef(ncid);
% %    netcdf.putVar(ncid,varid,sensor_matrix);
% %    netcdf.close(ncid);
% %    ncid2 = netcdf.open(filename,'NC_NOWRITE');
% %    data_copy = netcdf.getVar(ncid2,0);
% %    if isequal(sensor_matrix,data_copy)
% %        disp('Data match');
% %    else
% %        disp('Data mis-match');
% %    end
% 
% 
%     sensor_matrix = display_sensors_matrix_MAD(mask,sens(1:250),'r');
%     %imagesc(sensor_matrix)
% 
%     numrow = 101;
%     numcol = 143;
%     filename = ['figures/final_20220920_race_STL/netcdf/netcdf_top250_' num2str(Gamma) '.nc'];
%     ncid = netcdf.create(filename,'NC_WRITE');
%     dimidrow = netcdf.defDim(ncid,'rows',numrow);
%     dimidcol = netcdf.defDim(ncid,'length',numcol);
%     varid = netcdf.defVar(ncid,'loc_STL','NC_DOUBLE',[dimidrow dimidcol]);
%     netcdf.endDef(ncid);
%     netcdf.putVar(ncid,varid,sensor_matrix);
%     netcdf.close(ncid);
%     ncid2 = netcdf.open(filename,'NC_NOWRITE');
%     data_copy = netcdf.getVar(ncid2,0);
%     if isequal(sensor_matrix,data_copy)
%         disp('Data match');
%     else
%         disp('Data mis-match');
%     end
% 
% 
% 
%     sensor_matrix = display_sensors_matrix_MAD(mask,sens(1:200),'r');
%     %imagesc(sensor_matrix)
% 
%     numrow = 101;
%     numcol = 143;
%     filename = ['figures/final_20220920_race_STL/netcdf/netcdf_top200_' num2str(Gamma) '.nc'];
%     ncid = netcdf.create(filename,'NC_WRITE');
%     dimidrow = netcdf.defDim(ncid,'rows',numrow);
%     dimidcol = netcdf.defDim(ncid,'length',numcol);
%     varid = netcdf.defVar(ncid,'loc_STL','NC_DOUBLE',[dimidrow dimidcol]);
%     netcdf.endDef(ncid);
%     netcdf.putVar(ncid,varid,sensor_matrix);
%     netcdf.close(ncid);
%     ncid2 = netcdf.open(filename,'NC_NOWRITE');
%     data_copy = netcdf.getVar(ncid2,0);
%     if isequal(sensor_matrix,data_copy)
%         disp('Data match');
%     else
%         disp('Data mis-match');
%     end
% 
% 
%     sensor_matrix = display_sensors_matrix_MAD(mask,sens(1:150),'r');
%     %imagesc(sensor_matrix)
% 
%     numrow = 101;
%     numcol = 143;
%     filename = ['figures/final_20220920_race_STL/netcdf/netcdf_top150_' num2str(Gamma) '.nc'];
%     ncid = netcdf.create(filename,'NC_WRITE');
%     dimidrow = netcdf.defDim(ncid,'rows',numrow);
%     dimidcol = netcdf.defDim(ncid,'length',numcol);
%     varid = netcdf.defVar(ncid,'loc_STL','NC_DOUBLE',[dimidrow dimidcol]);
%     netcdf.endDef(ncid);
%     netcdf.putVar(ncid,varid,sensor_matrix);
%     netcdf.close(ncid);
%     ncid2 = netcdf.open(filename,'NC_NOWRITE');
%     data_copy = netcdf.getVar(ncid2,0);
%     if isequal(sensor_matrix,data_copy)
%         disp('Data match');
%     else
%         disp('Data mis-match');
%     end
% 
% 
%     sensor_matrix = display_sensors_matrix_MAD(mask,sens(1:100),'r');
%     %imagesc(sensor_matrix)
% 
%     numrow = 101;
%     numcol = 143;
%     filename = ['figures/final_20220920_race_STL/netcdf/netcdf_top100_' num2str(Gamma) '.nc'];
%     ncid = netcdf.create(filename,'NC_WRITE');
%     dimidrow = netcdf.defDim(ncid,'rows',numrow);
%     dimidcol = netcdf.defDim(ncid,'length',numcol);
%     varid = netcdf.defVar(ncid,'loc_STL','NC_DOUBLE',[dimidrow dimidcol]);
%     netcdf.endDef(ncid);
%     netcdf.putVar(ncid,varid,sensor_matrix);
%     netcdf.close(ncid);
%     ncid2 = netcdf.open(filename,'NC_NOWRITE');
%     data_copy = netcdf.getVar(ncid2,0);
%     if isequal(sensor_matrix,data_copy)
%         disp('Data match');
%     else
%         disp('Data mis-match');
%     end



end


% %% display reconstructions of sea surface temperature from multiscale sensors
% % this section can only run one gamma value at a time
% 
% %mkelp: find best gamma value then run reconstructions!
% 
% [~,~,time,mask,~] = read_data_enso([datpath 'sst.wkmean.1990-present_2.nc'],...
%     [datpath 'lsmask.nc']);
% times = datetime(1800,1,1,0,0,0) + days(time(indt));
% rmse_errors = [];
% mpe_errors = [];
% 
% for i=1:size(Xrecon,2)
%     figure; 
%     t = times(i);
%     display_fig_sst(Xrecon(:,i),mask,[],[-1.8 35.6]);
%     [rmse,mpe] = rmse_mpe(Xrecon(:,i),dat(:,:,indt(i)),mask);
%     rmse_errors = [rmse_errors; rmse];
%     mpe_errors = [mpe_errors; mpe];
%     if PRINT_FIG
%         export_fig([figpath 'recon' num2str([month(t),year(t)],'%d_%d')],'-pdf');
%     end
%     title(num2str([month(t),year(t)],'%d/%d'));
% end

%% display sensors
% this section also only runs for one gamma value at a time
% 
%  figure;
%  
%  display_sensors_MAD(mask,sens,'r');
%  if PRINT_FIG
%      %export_fig([figpath 'recon_sensors0.pdf']);
%      print('figures/final_20220920_race_STL/recon_sensors_lastone_TEST', '-dpdf', '-fillpage');
% 
%  end


%% Constructs Cost and Error vs Gamma Plot, similar to one in Fig. 7 of sensors-cost-paper
% this section is appropriate only when we have data for multiple gamma values

% figure;
% yyaxis right;
% hold on;
% plot(Gammas,Errors,'r-','LineWidth',1.5) % Plots error values in red
% h = gca; h.YColor = 'r';
% alpha(0.5)
% ylabel('Error','Color','r','FontName','Times')
% 
% yyaxis left
% hold on
% plot(Gammas,Costs,'k-','LineWidth',1.5) % Plots gamma values in black
% h = gca; h.YColor = 'k';
% alpha(0.5)
% set(gca,'FontSize',12)
% xlabel('\gamma','FontName','Times')
% ylabel('Cost','FontName','Times')
% hold on
% 
% print('figures/final2_STL/sensor_COSTs', '-dpdf', '-fillpage');