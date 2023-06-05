%--------------------------------------------------------------------------
%
% Load the SMA mesh (with EIT reconstruction parameters)
%
% 
% 
%--------------------------------------------------------------------------
clear;clc;close all;
% Set parameters:
params.bg         = 0.2;    % background conductivity of reference saline (S/m)
params.SMA_type   = 3;      % 1 for Allaire TREIT, 2 for Sophie oshpark sciospec, 3 for ZIF, 4 for JC  cherry rigid board
params.loc        = 2;      % 2 - no center electrode (v13), 3 is no corner electrode ()
params.eitType    = 'real'; % EIT recon type: real, imag, mag, phase
params.tik        = 1e5;    % tikhonov factor
params.res        = [30 30 1];
params.dbg_flg    = 0;      % debug flag
params.burst      = 1;      % number of files per sample, use 1 for now - change to 3 and average later
params.Ne         = 32;     % number of electrodes
params.com_clrbr  = 0.4;      %0.4; % to set a common colorbar choose the bounds here, otherwise set to 0 *might be overwritten by data sets that have clrbars defined
params.iter       = 5;      % number of iterations for EIT recons
params.filter     = 1;      % filter = 1 to run Sophie's filter, Shannon had something with filter = 3?
params.dbg_flg    = 0;
params.nmFreqs    = 41;     % number of frequencies: 100 Hz to 1 MHz, 41 log-spaced frequencies
params.f_sel      = [1,11,21,28,31,41]; % select index of frequency, set to vector to loop through frequencies
params.f_test     = [100,1e3,1e4,50119,1e5,1e6]; % corresponding frequencies for f_sel
% params.f_test     = [100,1e3,1e4,5e4,1e5,1e6]; % corresponding frequencies for f_sel
params.iz_sel     = [1:4]; % select Z layer of recon, 5 total, 5 is the closest but
                  % has a lot of impact from electrode artifact -> don't use 
                  % for EIT images. Use 4 as the closest, 1 as the furthest
                  % set to vector to loop through layers
params.showRefVcomp = 1; % set 1 to show comparision of voltages between ref and test, set 0 otherwise
params.crp          = 1; % if = 1, crop corse_grid display of binary map/EIT recon
params.binMaps      = 1; % if 1 add binary map to Latex slide, if 0 ignore
params.sim          = 0; % if 1 then EIT recons for sim, otherwise experimental data

% ----------------- *** Select current injection *** ---------------------
% current loaded in with file set, uncomment if needed
% I = 1e-4; % 100 uA
I = 1e-3; % 1mA 
fprintf('Injected current: %f mA\n',I*1000)

%--------------------------------------------------------------------------
% Add paths
% addpath(genpath('S:/digihisto/Ethan/NDRM'))
% addpath('S:/digihisto/Ethan/EKM_utility')
% % addpath('S:/eit/Kossmann/Sciospec_EIT/mfiles')
% addpath('S:/eit/Allaire/Matlab_Functions')

addpath(genpath('/Volumes/jumbo/digihisto/Ethan/NDRM'))
addpath('/Volumes/jumbo/digihisto/Ethan/EKM_utility')
% addpath('S:/eit/Kossmann/Sciospec_EIT/mfiles')
addpath('/Volumes/jumbo/eit/Allaire/Matlab_Functions')
addpath('/Volumes/jumbo/digihisto/Ethan/EKM_utility/meshing')
addpath(genpath('/Volumes/jumbo/eit/Allaire/Zenia_meat_database'))
%%
%--------------------------------------------------------------------------
% Load the mesh
%load dat/Green2D_1p5
disp('loading mesh...')

helec    = 0.1;           % // h-value on electrodes
hfine    = 0.3;          % // h-value in the fine area
hfar     = 5;
iter = 5;
geo_str  = ['probe_rounded2x2Ipairs_larger_25_5x5Vs','_h',ifdec(num2str(helec)), ...
    '_',ifdec(num2str(hfine)),'_',ifdec(num2str(hfar))];
eval(['load optmesh',num2str(iter),'_',geo_str,'_large_circleEs_sk msh'])
%%
%--------------------------------------------------------------------------
% plot example image
% load data
pathfront = '/Volumes/jumbo/eit/Allaire/Zenia_meat_database/';
data = 5;
typeset = {'Muscle','Fat','Fat Strip'};
locset = {'Loc1','Loc5','Loc6'};

switch(data)
   case 1 
      load([pathfront '2023_03_24_muscle_loc1_22g/muscle_loc1_22g_real_cropped_100_Hz_z4.mat'])
      type = typeset{1};
      loc = locset{1};
      lim = [-0.3,0.3];
      filename = 'muscle_loc1_22g_real_cropped_100_Hz_z';
   case 2 
       load([pathfront '2023_03_24_fat_loc1_21g/fat_loc1_21g_real_cropped_100_Hz_z4.mat'])
       type = typeset{2};
      loc = locset{1};
      lim = [-1.5,1.5];
      filename = 'fat_loc1_21g_real_cropped_100_Hz_z';
   case 3
      load([pathfront '2023_03_24_fatstrip_loc1_25g/fat_strip_loc1_25g_real_cropped_100_Hz_z4.mat'])
      type = typeset{3};
      loc = locset{1};
      lim = [-0.3,0.3];
      filename = 'fat_strip_loc1_25g_real_cropped_100_Hz_z';
   case 4
      load([pathfront '2023_03_24_fatstrip_loc5c_29g/fat_strip_loc5c_29g_real_cropped_100_Hz_z4.mat'])
      type = typeset{3};
      loc = locset{2};
      lim = [-0.4,0.4];
      filename = 'fat_strip_loc5c_29g_real_cropped_100_Hz_z';
   case 5 
      load([pathfront '2023_03_24_fatstrip_loc6c_27g/fat_strip_loc6c_27g_real_cropped_100_Hz_z4.mat'])
      type = typeset{3};
      loc = locset{3};
      lim = [-0.4,0.4];
      filename = 'fat_strip_loc6c_27g_real_cropped_100_Hz_z';
   otherwise
      fprintf('Invalid data\n' );
end

disp('test')

for j= 1:4 % loop through Z layers
    if params.crp == 1
    %   disp('Crop corse_grid');
    % ignore line below, corse_grid cropped when data was saved - you are
    % loading the output of this function
    % [pltset_xy_crp] = cropCorseGrid(pltset_xy); %crop corse_grid
        rectsx = pltset_xy_crp(j).rectsx;
        rectsy = pltset_xy_crp(j).rectsy;
        inds   = pltset_xy_crp(j).inds;
    else
        rectsx = pltset_xy(j).rectsx;
        rectsy = pltset_xy(j).rectsy;
        inds   = pltset_xy(j).inds;
    end
    figure; hold on;
    %patch(rectsx,rectsy,de_s,'linestyle','none')
    if strcmp(params.eitType,'real')
        sigmas_idx = sigmas(inds);
        patch(rectsx,rectsy,0*rectsx,sigmas_idx,'linestyle','none')% ,'facealpha',0.8)
    %             patch(pltset_xy(j).rectsx, ... 
    %                     pltset_xy(j).rectsy, ...
    %                     0*pltset_xy(j).rectsx, ...
    %                     de_s(pltset_xy(j).inds),'linestyle','none')% ,'facealpha',0.8)
    elseif strcmp(params.eitType,'imag')
        epss_idx   = epss(inds);
        patch(rectsx,rectsy,0*rectsx,epss_idx,'linestyle','none')% ,'facealpha',0.8)
    end
    view(2)
    disp('test')

    %---------------------------------------------------------------------------
    % Plot the electrodes
    %plot_msh_elecs_only(msh,1000,1,1)
    i_elecs = [1:33];
    ebnds = plot_msh_elecs_outline(msh,1000,i_elecs,0); % last variable input 0 to not display electrode number, 1 to display electrode number
    view(2)
    % plot max sensitivity colorbar for each image - can change to be
    % common colorbar across Z and frequency
%     if strcmp(params.eitType,'real') 
%         caxis(max(abs(sigmas_idx))*[-1 1])
%     elseif strcmp(params.eitType,'imag')
%         caxis(max(abs(epss_idx))*[-1 1])
%     end

    caxis(lim)
    figFolder = 'figs';
    
    axis off
    axis equal 
    colorbar;
    inc_trace    = 0; %1 if you want to trace an outline of the inclusion
    inc_cent     = [0,0];
    inc_diam     = 2; %in mm
    if inc_trace == 1
        plot(inc_cent(1)+inc_diam/2*cos(ts),inc_cent(2)+inc_diam/2*sin(ts),'--g')
    end
    % ff = round(ftest(i)); % frequency of recon
    % t = [num2str(ff/1000),' kHz, z = ',num2str(j)];
   
    t = [type,' ',loc, ' 100 Hz Z=',num2str(j),' Real Difference EIT Recon'];
    title(t)
    
    dest_dir = ['/Users/zenia/Desktop/ENGS 111/final_project/eit_images'];
    fname = [filename,num2str(j),'.png'];
    
    % saveas(fig,filename)
    fig_file = fullfile(dest_dir , fname);
    
    %saveas(a, filename)  %save the file there directory
    saveas(gcf, fig_file)

end