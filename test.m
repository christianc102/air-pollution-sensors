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
    display_sensors_MAD(mask,sens,'r');
    if PRINT_FIG
        filename = ['figures/final2_STL/recon_sensors_ideal_' num2str(Gamma)];
        %export_fig([figpath filename]);
        print(filename, '-dpdf', '-fillpage');

        %export_fig([figpath 'recon_sensors.pdf']);
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

    figure;
    display_sensors_MAD(mask,sens,'r');
    if PRINT_FIG
        filename = ['figures/final_20220920_race_STL/recon_sensors_top250_' num2str(Gamma)];
        print(filename, '-dpdf', '-fillpage');
    end


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

 figure;
 
 display_sensors_MAD(mask,sens,'r');
 if PRINT_FIG
     %export_fig([figpath 'recon_sensors0.pdf']);
     print('figures/final_20220920_race_STL/recon_sensors_lastone_TEST', '-dpdf', '-fillpage');

 end


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