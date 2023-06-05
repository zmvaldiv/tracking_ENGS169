%--------------------------------------------------------------------------
%
% Run a difference reconstruction on two sets of measured Sciospec data
%
% 
% 
%--------------------------------------------------------------------------
clear;clc;close all

%--------------------------------------------------------------------------
% Set parameters:
tik          = 1e5;                         % Tikhonov Reg. param values
SMA_type     = 3; %1 for Allaire TREIT, 2 for Sophie oshpark sciospec, 3 for ZIF, 4 for JC  cherry rigid board
loc          = 2; %2 - no center electrode (v13), 3 is no corner electrode ()
sel_case     = 'ENGS57Demo';
% datpth0      = 'S:/eit/Allaire/ex_vivo_probe/PR002_2023_01_25'; %*must use / not \
com_clrbr    = 0%0.4; % to set a common colorbar choose the bounds here, otherwise set to 0 *might be overwritten by data sets that have clrbars defined
        
switch sel_case
    case '2022_09_15_meat_1mA'
        datpth0      = 'S:/eit/Allaire/ex_vivo_probe/2022_09_15_meat'; %*must use / not \
        load('2022_09_15_meat/2022_09_15_meat_files.mat') 
    case 'PR002_100uA'
        datpth0      = 'S:/eit/Allaire/ex_vivo_probe/PR002_2023_01_25'; %*must use / not \
        load([datpth0,'/PR002_100uA_files.mat']) 
    case 'PR002_1mA'
        datpth0      = 'S:/eit/Allaire/ex_vivo_probe/PR002_2023_01_25'; %*must use / not \
        load([datpth0,'/PR002_1mA_files.mat']) 
    case '2023_03_24_meat_1mA'
        datpth0      = 'S:/eit/Allaire/ex_vivo_probe/2023_03_24_meat'; %*must use / not \
%         load([datpth0,'/exVivoBovine_allfiles.mat']) % all files

%         load([datpth0,'/exVivoBovine_regfiles.mat']) % unmodified files
%         com_clrbr    = [0.45,0.55,0.4,0.5,0.5,2.5,5,0.4,0.4,1.1,3.5,0.5,0.5,0.5,0.6,1.2,0.5,10,3,2.5,3.5,0.6,0.4,0.8,0.6];
%         I=reg_Is_1mA(1);
%         file_nms=reg_tstFiles_1mA;
%         ref_nm=reg_refFiles_1mA{1};
%         load([datpth0,'/bovine_fat_homo_1mA_files.mat']) % 100% fat files
%         load([datpth0,'/bovine_muscle_homo_1mA_files.mat']) % 100% muscle files
%         load([datpth0,'/bovine_muscle_fat_hetero_1mA_files.mat']) % muscle/fat combo files
% 
        load([datpth0,'/exVivoBovine_modfiles.mat']) % modified files
        com_clrbr    = [0.3,0.3,0.4,0.3,0.3,0.45,1,0.3,0.4,0.5,2.5,3,4,0.4,0.8,0.4,0.45];
        I=mod_Is_1mA(1);
        file_nms=mod_tstFiles_1mA;
        ref_nm=mod_refFiles_1mA{1};
%         load([datpth0,'/bovine_saline0o2_1mA_files.mat']) % 0.2S/m saline files
%         load([datpth0,'/bovine_saline0o7_1mA_files.mat']) % 0.7S/m saline files
%         load([datpth0,'/bovine_cautery_1mA_files.mat']) % cautery files
%         load([datpth0,'/bovine_dye_1mA_files.mat']) % methylene blue dye files
%         load([datpth0,'/bovine_lidocaine_1mA_files.mat']) % cautery files
    case 'ENGS57Demo'
         datpth0      = 'S:\eit\Zenia\engs169_test_20230604'; %*must use / not \
         I = 1e-4; % 1mA 
         %file_nms={'\20230604 11.38.08\setup_full_00002\setup_full_00002_3.eit'};
         file_nms={'\20230604 12.48.22\setup_full_00002\setup_full_00002_3.eit'};
         ref_nm='setup_full_00001_1.eit';
end

% using loop for file_nms below, uncomment to use single file
% test_nm       = 'muscle1_45g_hh\setup_00001\setup_00001_1.eit';
% ref_nm       = 'blank_0o2_100uA/setup_00002/setup_00002_1.eit'; %
% reference loaded in with file set, uncomment if needed

% ----------------- *** Select current injection *** ---------------------
% current loaded in with file set, uncomment if needed
% I = 1e-4; % 100 uA
% I = 1e-3; % 1mA 
fprintf('Injected current: %f mA\n',I*1000)
%--------------------------------------------------------------------------
rmv_elecs    = [];
c_max        = []; %colorbar limit if known, otherwise, empty set
mx           = 0;
real_flg     = 1;                           % Real flag
bg           = 0.2;                         % Background conductivity
enum         = 33;
outer        = 8;
iter         = 5;
inc_trace    = 0; %1 if you want to trace an outline of the inclusion
inc_cent     = [0,0];
inc_diam     = 2; %in mm
ts           = linspace(0,2*pi,100);
filter       = 3;
dbg_flg      = 0;
% f_sel        = [1,11,21,31,41]; % select index of frequency, set to vector to loop through frequencies
% f_sel        = [1,11,21,28,31,41]; % select index of frequency, set to vector to loop through frequencies
f_sel        = [1,11,21,28,31]; % select index of frequency, set to vector to loop through frequencies
iz_sel       = [1:4]; % select Z layer of recon, 5 total, 5 is the closest but
                  % has a lot of impact from electrode artifact -> don't use 
                  % for EIT images. Use 4 as the closest, 1 as the furthest
                  % set to vector to loop through layers
showRefVcomp = 1; % set 1 to show comparision of voltages between ref and test, set 0 otherwise
crp          = 1; % if = 1, crop corse_grid display of binary map/EIT recon

%--------------------------------------------------------------------------
% Add paths
addpath(genpath('S:/digihisto/Ethan/NDRM'))
addpath('S:/digihisto/Ethan/EKM_utility')
% addpath('S:/eit/Kossmann/Sciospec_EIT/mfiles')
addpath('S:/eit/Allaire/Matlab_Functions')
addpath('S:\eit\Allaire\ex_vivo_probe\dat')

%--------------------------------------------------------------------------
% Load the mesh
%load dat/Green2D_1p5
disp('loading mesh...')
if outer == 8
    helec    = 0.1;           % // h-value on electrodes
    hfine    = 0.3;          % // h-value in the fine area
    hfar     = 5;
    geo_str  = ['probe_rounded2x2Ipairs_larger_25_5x5Vs','_h',ifdec(num2str(helec)), ...
        '_',ifdec(num2str(hfine)),'_',ifdec(num2str(hfar))];
    eval(['load optmesh5_probe_rounded2x2Ipairs_larger_25_5x5Vs_h0o1_0o3_5_large_circleEs_sk.mat'])
end
%--------------------------------------------------------------------------
% Load the Jacobian data
disp('loading Jacobian...')
% if SMA_type == 1
%     eval(['load dat/J_8outerallI_sens1ovr50_bg',ifdec(num2str(bg)), ...
%         '_res30x30x20_h0o2_0o4_5_opt5_Sciospeciivv_loc',num2str(loc),' J L corse2fine corse_grid Zref'])
% elseif SMA_type == 2
%     if loc == 3
%         eval(['load dat/J_8outerallI_sens1ovr50_bg',ifdec(num2str(bg)),...
%             '_res30x30x20_h0o2_0o4_5_opt5_Sciospeciivv_Sophie_loc',num2str(loc)])
%     elseif loc == 2
%         load dat/J_8outerallI_sens1ovr50_bg0o2_res30x30x20_h0o2_0o4_5_opt5_Sciospeciivv_ZIF_loc3 % ZIF doesn't use center, so this should match
%     end
% elseif SMA_type == 3
%         load dat/J_8outerallI_sens1ovr50_bg0o2_res30x30x20_h0o2_0o4_5_opt5_Sciospeciivv_ZIF_loc3
% elseif SMA_type == 4
%         load dat/J_8outerallI_sens1ovr50_bg0o2_res30x30x20_h0o2_0o4_5_opt5_Sciospeciivv_ZIF_loc3
%         %Jacobian is in terms of electrode number without V13, no changes needed for JC cherry board
% end 

load J_8outerallI_sens1ovr50_bg0o2_res30x30x20_h0o2_0o4_5_opt5_Sciospeciivv_ZIF_loc3.mat

for n = 1:size(file_nms,1) % loop through files
    test_nm = file_nms{n};
    disp(['processing: ',test_nm])
    %--------------------------------------------------------------------------
    % Read the Sciospec data
    if SMA_type == 1 % Allaire TREIT
        % IIVV patterns
        iivv = sciospec_eit_iivv(loc); iivv0 = iivv;
        iis = iis_sciospec(); iis0 = iis;
    elseif SMA_type == 2 %Sophie oshpark sciospec
        % Read the Sciospec data
        [iis0,vs_re_test0,vs_im_test0,ftest] = read_sciospecEIT_dat(datpth0,test_nm);
        [iis0,vs_re_ref0,vs_im_ref0,fref] = read_sciospecEIT_dat(datpth0,ref_nm);
        vs_re_test = scio2mesh_sophie(vs_re_test0,loc);
        vs_im_test = scio2mesh_sophie(vs_im_test0,loc);
        vs_re_ref = scio2mesh_sophie(vs_re_ref0,loc);
        vs_im_ref = scio2mesh_sophie(vs_im_ref0,loc);
        % IIVV patterns
        load dat/ZIFskip0iivvr13.mat % SMA 8 current drive, 24 voltage pick up (no electrode v13)
                                     % loads iivv patterns (7728 x 4 double)
        if loc == 2
            load dat/sophie_oshpark_corner_map_plg2plg.mat
        %         map_plg2plg = sciospec_Sophie_SMA_map(); % old function used to create mapping
        elseif loc == 3
            load dat/sophie_oshpark_center_map_plg2plg.mat
        end
        % old code
        %     iivv0 = sciospec_num_32_iivv();
        %     iis0 = iis; %unique(iivv0(:,1:2),'rows','stable');
        %     map_plg2plg = sciospec_Sophie_SMA_map();
        %     for n = 1:4
        %         iivv(:,n) = map_plg2plg(iivv0(:,n),loc);
        %     end
        % %     iis = unique(iivv(:,1:2),'rows','stable');
    elseif SMA_type == 3 % ZIF electrode array with small connectors
        % Read the Sciospec data
        [iis0,vs_re_test0,vs_im_test0,ftest] = read_sciospecEIT_dat(datpth0,test_nm);
        [iis0,vs_re_ref0,vs_im_ref0,fref] = read_sciospecEIT_dat(datpth0,ref_nm);
        vs_re_test = scio2mesh_ZIF(vs_re_test0);
        vs_im_test = scio2mesh_ZIF(vs_im_test0);
        vs_re_ref = scio2mesh_ZIF(vs_re_ref0);
        vs_im_ref = scio2mesh_ZIF(vs_im_ref0);
        % IIVV patterns
        load dat/ZIFskip0iivvr13.mat % SMA 8 current drive, 24 voltage pick up (no electrode v13)
                                     % loads iivv patterns (7728 x 4 double)
        load dat/ZIFmap_plg2plg.mat  % column 1 is sciospec channel, column 2 is the corresponding electrode
    elseif SMA_type == 4 %JC Cherry
        % Read the Sciospec data
        [iis0,vs_re_test0,vs_im_test0,ftest] = read_sciospecEIT_dat(datpth0,test_nm);
        [iis0,vs_re_ref0,vs_im_ref0,fref] = read_sciospecEIT_dat(datpth0,ref_nm);
        vs_re_test = scio2mesh_JCcherry(vs_re_test0);
        vs_im_test = scio2mesh_JCcherry(vs_im_test0);
        vs_re_ref = scio2mesh_JCcherry(vs_re_ref0);
        vs_im_ref = scio2mesh_JCcherry(vs_im_ref0);
        % IIVV patterns -     %%% NO EDITS MADE YET
        load dat/ZIFskip0iivvr13.mat % SMA 8 current drive, 24 voltage pick up (no electrode v13)
                                     % loads iivv patterns (7728 x 4 double)
                                     % exhaustive iivv patterns in terms of electrode number without V13, 
                                     % no changes needed for JC cherry board
        load dat/JCcherrymap_plg2plg.mat  % column 1 is sciospec channel, column 2 is the corresponding electrode
    else
    %     [iis,vs_re_test,vs_im_test,ftest] = read_sciospecEIT_dat(datpth0,test_nm);
    %     [iis,vs_re_ref,vs_im_ref,fref] = read_sciospecEIT_dat(datpth0,ref_nm);
    end
    
    %--------------------------------------------------------------------------
    % IIVV patterns
    %     iis(:,1) = map_plg2plg(iis0(:,1),loc);
    %     iis(:,2) = map_plg2plg(iis0(:,2),loc);
        iis(:,1) = map_plg2plg(iis0(:,1),2); % col 1 is I+
        iis(:,2) = map_plg2plg(iis0(:,2),2); % col 2 is I-
    
    c2f = corse2fine;
    
    %--------------------------------------------------------------------------
    for i = f_sel % loop through selected frequencies
        for j = iz_sel % loop through selected layers
            %convert Sciospec data into Z data (first freq)
            vs_test = vs_re_test(:,:,i); vs_ref = vs_re_ref(:,:,i);
            if filter == 2
               [c_ind,e_ind] = find(abs(vs_test) < 0.001); 
               [c_ind2,e_ind2] = find(abs(vs_ref) < 0.001);
               rmv_elecs = [e_ind;e_ind2];
               rmv_elecs = unique(rmv_elecs);
            end
            urmv_pats = [];
            Zs_test = calcIIVV_from_single_ended_iis(iis,iivv,vs_test/I);
            Zs_ref  = calcIIVV_from_single_ended_iis(iis,iivv,vs_ref/I);
    %         if SMA_type == 2
    %             Zs_test  = calcIIVV_from_single_ended_iis(iis,iivv,vs_test/I);
    %             Zs_ref  = calcIIVV_from_single_ended_iis(iis,iivv,vs_ref/I);
    % %             Zs_test  = calcIIVV_from_single_ended_iis(iis0,iivv0,vs_test/I);
    % %             Zs_ref  = calcIIVV_from_single_ended_iis(iis0,iivv0,vs_ref/I);
    %         end
    
            %Filter by removing Zs with |Z|<0.01
            if filter == 1
               rmv_pats = find(abs(Zs_test./Zs_ref)>10);
               rmv_pats_low = find(abs(Zs_test./Zs_ref < 0.1));
               urmv_pats = [rmv_pats_low;rmv_pats]; urmv_pats = unique(urmv_pats);
            %    rmv_test = find(abs(Zs_test) < 0.01); %Zs_meas(indx) = NaN; 
            %    rmv_ref = find(abs(Zs_ref) < 0.01);
            %    urmv_pats = [rmv_test;rmv_ref]; urmv_pats = unique(urmv_pats);
            end
            if filter ~= 0 
                if ~isempty(rmv_elecs)
                    rmv_list = [];
                    for a = 1:length(rmv_elecs)
                        [indx,b] = find(iivv == rmv_elecs(a));
                        rmv_list = [rmv_list; indx];
                    end
                    urmv_pats = [urmv_pats;rmv_list];
                    urmv_pats = unique(urmv_pats);
                end
                Zs_test(urmv_pats,:) = [];
                Zs_ref(urmv_pats,:)  = [];
                J(urmv_pats,:)      = [];
                Ztestav   = Zs_test;
                Zrefav    = Zs_ref;
    
            %     n=1;
            %     for i=1:length(iivv)
            %         if ~ismember(i,rmv)
            %             Ztestav(n) = Zs_test(i);
            %             Zrefav(n) = Zs_ref(i);
            %             iivv_s(n,:) = iivv(i,:);
            %             n = n+1;
            %         end
            %     end
            else
                Ztestav   = Zs_test;
                Zrefav    = Zs_ref;
            end
    
            %--------------------------------------------------------------------------
            % Construct the patch matrices
            uxs = unique(corse_grid.node(1:end-1,1));
            uys = unique(corse_grid.node(1:end-1,2));
            uzs = unique(corse_grid.node(1:end-1,3));
            [corse_grid.xs,ix] = intersect(corse_grid.xs,uxs);
            [corse_grid.ys,iy] = intersect(corse_grid.ys,uys);
            [corse_grid.zs,iz] = intersect(corse_grid.zs,uzs);
            corse_grid.dxs     = corse_grid.dxs(ix);
            corse_grid.dys     = corse_grid.dys(iy);
            corse_grid.dzs     = corse_grid.dzs(iz);
            % mesh units in meters
    
            [pltset_xy,uxs,uys,uzs,pltset_xz,pltset_yz] = constr_corse_grid_patches_varydxyzs(corse_grid,1000,1);
            if crp == 1
%                 disp('Crop corse_grid');
                [pltset_xy_crp] = cropCorseGrid(pltset_xy); %crop corse_grid
                rectsx = pltset_xy_crp(j).rectsx;
                rectsy = pltset_xy_crp(j).rectsy;
                inds   = pltset_xy_crp(j).inds;
            else
                rectsx = pltset_xy(j).rectsx;
                rectsy = pltset_xy(j).rectsy;
                inds   = pltset_xy(j).inds;
            end
    
            %--------------------------------------------------------------------------
            % Calculate the difference in impedances
            Z_diff = Ztestav - Zrefav;
            % Absolute conductivities
            %Z_diff = Zrefav;
            %---------------------------------------------------------------------------
            % Calculate the reconstruction
            % Z_diff = Z_diff';
            de_s         = ( (J'*J) + tik * (L'*L) ) \ (J'*Z_diff);
    
            Q     = tik*((L')*L);
            conds = (J'*J+Q)\(J'*Z_diff );
            conds = real(conds); %and this also real part only (conds)
            de_s  = conds;
            de_s_idx = conds(inds);
            
            figure
            hold on
            %patch(rectsx,rectsy,de_s,'linestyle','none')
            patch(rectsx,rectsy,0*rectsx,de_s_idx,'linestyle','none')% ,'facealpha',0.8)
    %             patch(pltset_xy(j).rectsx, ... 
    %                     pltset_xy(j).rectsy, ...
    %                     0*pltset_xy(j).rectsx, ...
    %                     conds(pltset_xy(j).inds),'linestyle','none')% ,'facealpha',0.8)
            view(2)
            %---------------------------------------------------------------------------
            % Plot the electrodes
            %plot_msh_elecs_only(msh,1000,1,1)
            i_elecs = [1:33];
            ebnds = plot_msh_elecs_outline(msh,1000,i_elecs,0); % last variable input 0 to not display electrode number, 1 to display electrode number
            view(2)
            if length(com_clrbr) == 1 && ~sum(com_clrbr) == 0
                caxis(com_clrbr*[-1 1]) %one common value
            elseif sum(com_clrbr) == 0
                caxis(max(abs(de_s_idx))*[-1 1])
            else
                caxis(com_clrbr(n)*[-1 1]) %index value from vector
            end
            %caxis([-.06 .1])
            axis off
            axis equal 
            colorbar
    
            if inc_trace == 1
                plot(inc_cent(1)+inc_diam/2*cos(ts),inc_cent(2)+inc_diam/2*sin(ts),'--g')
            end
            ff = round(ftest(i)); % frequency of recon
            t = [num2str(ff/1000),' kHz, z = ',num2str(j)];
            title(t)
            if crp == 1
                % save figure
                saveas(gcf,[datpth0,'\figs',test_nm(1:end-30),'_cropped_',num2str(ff),'_','Hz','_z',num2str(j)],'.png')
                % save EIT data
                eval(['save ',datpth0,'/dat/',test_nm(1:end-30),'_cropped_',num2str(ff),...
                      '_','Hz','_z',num2str(j),' pltset_xy_crp conds'])   
            else
                % save figure
                saveas(gcf,[datpth0,'\figs',test_nm(1:end-30),'_',num2str(ff),'_','Hz','_z',num2str(j)],'.png')
                % save EIT data
                eval(['save ',datpth0,'\dat',test_nm(1:end-30),'_',num2str(ff),...
                      '_','Hz','_z',num2str(j),' pltset_xy conds'])   
            end
            
            
        end
        if showRefVcomp ==1
            figure; hold on;
            vs_re_T = vs_test';
            V_ref_T = vs_ref';
            
            plot(vs_re_T(:))
            plot(V_ref_T(:))
    
            lbl_fmt_fig('Concatenated Electrode Data','Voltage (V)',[num2str(round(ff/1000)),'_','kHz'],'','',12)
            legend('V_{incl}','V_{blank}')
            hold off
            % set(gca,'xlim',[1 32])
            %     saveas(gcf,['figs/example_meas_sim_comp_all_iipats_default_type'],'png')
            
            % save figure
            saveas(gcf,[datpth0,'\figs',test_nm(1:end-30),'_',num2str(ff),'_','Hz','_concat_patterns'],'.png')
    
            figure
            hold on
            plot(Zs_test(:),'.k','markersize',8);
            plot(Zs_ref,'.r','markersize',6);
            legend('Z{Meas}','Z{Ref}')
            lbl_fmt_fig('IIVV','Impedance (\Omega)','Z_{tissue} vs Z_{blank}','','',12)
            % save figure
            saveas(gcf,[datpth0,'\figs',test_nm(1:end-30),'_',num2str(ff),'_','Hz','_Zcomp'],'.png')
        end
    end
    close all; %close all figures to prevent matlab from crashing
end