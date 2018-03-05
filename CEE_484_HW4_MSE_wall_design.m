clear all
close all
clc

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Introduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Author: Enok Cheon, Ph.D. Student in UIUC (Supervisor: Prof Timothy D. Stark)
Date: Completed on 28th Feb, 2018
Description:
    The script was written to primarily solve CEE484 HW4 problem, which is design of MSE walls.
    The script calculates minimum length of geogrid required to resist loads from the pressure from 
    (1) earth and (2) surcharge with groundwater (GW) level is low to not effect the MSE wall.
Assumption and Considerations:
    (1) groundwater (GW) level is very low - does not effect the MSE wall stability.
    (2) global stability is adequate.
    (3) the surcharge load is uniform surcharge; no point load or line load
    (4) geogrid is assumed to be extandable; hence, follows failure wedge defined by Rankine theory
        for internal stability check. 
    (5) all the geogrid length and vertical spacing is kept constant for simplicity.
    (6) pull-out resistance factor of geogrid is assumed to be 0.7*tan(effective friction angle of MSE fill)
    (7) inclination of wall and backfill is vertical and horizontal, respectively.
    (8) minimum width-to-height ratio (B/H) for MSE wall is 0.7
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wall 
H = 24;     % wall height - unit: ft
Bmin = 0.7*H;
Bmax = 24;
dB = 0.25;

%% soil properties
% fill for MSE
gamma_f = 115;	% unit weight - unit: pcf
phi_f = 32;	    % effective friction angle
c_f = 0;		% effective cohesion - unit: psf

% soil behind wall
gamma_bw = 113;  % unit weight - unit: pcf
phi_bw = 30;     % effective friction angle
c_bw = 0;        % effective cohesion - unit: psf

% soil underneath wall
q_ult_uw = 12000;  % ultimate bearing capacity - unit: pcf

%% unform surcharge
q = 300;	% unit: psf

%% geogrid
Cr = 0.8;	% coverage ratio 
Tult = 12500; % ultimate tensile strength - unit: psf
RF = 12.5/3; % reduction factor for tensiel strength
F_pull = 0.7*tand(phi_f); % pullout resistance factor of geogrid
scale_effect = 0.8; 
spacing_v_min = 2;  % vertical spacing of geogrids - unit: ft
spacing_v_max = 4;  % vertical spacing of geogrids - unit: ft
spacing_v_d = 0.5;  % vertical spacing of geogrids - unit: ft

%% FS
FS_o_lim = 2;		% FS requirement for overturning
FS_bc_lim = 2.5;	% FS requirement for bearing capacity
FS_s_lim = 1.5;		% FS requirement for sliding
FS_r_lim = 1.2;		% FS requirement for rupture
FS_p_lim = 1.5;		% FS requirement for pullout

%% format of numbering
format short
significantFigure = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allsummaryData = [];
indiv_geogrid_Dat = {};

b = [Bmin, ceil(Bmin):dB:Bmax];
sv = [spacing_v_min:spacing_v_d:spacing_v_max];

rowN = 1;
rowNN = 1;
for i = 1:length(b)
    for ii = 1:length(sv) 
	
    	B = b(i); 	    % width
        SV = sv(ii);     % vertical spacing

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%% external stability
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%% force and moment
    	% vertical force
    	W_fill = gamma_f*B*H;	
    	Q = q*B;
    	sum_Fv = W_fill + Q;

    	% horizontal force
    	Ka = (1-sind(phi_bw))/(1+sind(phi_bw));
    	Pa_s = 0.5*gamma_bw*Ka*H^2;
    	Pa_q = q*H*Ka;
    	sum_Fh = Pa_s + Pa_q;

    	%% overturning analysis
    	Mo = Pa_q*(H/2) + Pa_s*(H/3);
    	Mr = (W_fill+Q)*(B/2);
    	
    	FS_o = Mr/Mo;
    	if FS_o >= FS_o_lim
    		check_o = 1;
    	else 
    		check_o = 0;
    	end

    	%% bearing capcity
    	% resultant moment about the MSE wall heel
    	sum_M = sum_Fv*(B/2) + Mo;
    	e = (sum_M/sum_Fv)-(B/2);
    	if e <= B/6
    		% q_max = (sum_Fv/B)*(1+(6*e/B));
            q_max = sum_Fv/(B-2*e);
    	else 
    		q_max = 99999999;
    	end

    	FS_bc = q_ult_uw/q_max;
    	if FS_bc >= FS_bc_lim
    		check_bc = 1;
    	else 
    		check_bc = 0;
    	end

    	%% sliding
    	interface_delta = min([phi_f, phi_bw]);
    	Fh_r = sum_Fv*tand(interface_delta);

    	FS_s = Fh_r/sum_Fh;
    	if FS_s >= FS_s_lim
    		check_s = 1;
    	else 
    		check_s = 0;
    	end

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%% internal stability
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%% Tmax and resisting zone length for each geogrid elevation
    	indiv_geogrid = [];
    	Z_num = ceil(H/SV);
        startingZ = H-(Z_num-1)*SV;
    	indiv_geogrid(:,1) = [startingZ:SV:H];
    	
    	% Tmax
    	Kr = (1-sind(phi_f))/(1+sind(phi_f)); 	% Ka = Kr for geogrids
    	indiv_geogrid(:,2) = indiv_geogrid(:,1)*Kr*gamma_f;	
    	indiv_geogrid(:,3) = indiv_geogrid(:,2)*SV/Cr;	% Tmax

    	% resisting zone length
    	indiv_geogrid(:,4) = (H*ones(Z_num,1) - indiv_geogrid(:,1))*tand(45-(phi_f/2)); % Lh
    	indiv_geogrid(:,5) = H*ones(Z_num,1) - indiv_geogrid(:,4); %Lr

    	%% rupture- long term
    	indiv_geogrid(:,6) = (Tult/RF)*ones(Z_num,1);
    	for loop_rup = 1:Z_num
    		indiv_geogrid(loop_rup,7) = indiv_geogrid(loop_rup,6)/indiv_geogrid(loop_rup,3);
    		if indiv_geogrid(loop_rup,7) >= FS_r_lim
    			indiv_geogrid(loop_rup,8) = 1;
    		else 
    			indiv_geogrid(loop_rup,8) = 0;
    		end
    	end

    	%% pullout 
    	for loop_pull = 1:Z_num
            indiv_geogrid(loop_pull,9) = 2*(indiv_geogrid(loop_pull,1)*gamma_f)*F_pull*Cr*scale_effect*(indiv_geogrid(loop_pull,5)-1);
    		indiv_geogrid(loop_pull,10) = indiv_geogrid(loop_pull,9)/indiv_geogrid(loop_pull,3);
    		if indiv_geogrid(loop_pull,10) >= FS_p_lim
    			indiv_geogrid(loop_pull,11) = 1;
    		else 
    			indiv_geogrid(loop_pull,11) = 0;
    		end
    	end

    	%% overall internal stability FS check
    	% rupture
    	if sum(indiv_geogrid(:,8)) == Z_num
    		check_r = 1;
    	elseif sum(indiv_geogrid(:,8)) < Z_num 
    		check_r = 0;
    	end

    	% pullout
    	if sum(indiv_geogrid(:,11)) == Z_num
    		check_p = 1;
    	elseif sum(indiv_geogrid(:,11)) < Z_num
    		check_p = 0;
    	end

    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%% summary
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%% overall stability FS check
    	if check_o + check_bc + check_s + check_r + check_p == 5 
    		overall_check = 1;
    	else 
    		overall_check = 0;
    	end

        %% amount of geogrids used
        overallCost = Z_num*B;

    	%% collect data
    	allsummaryData(rowN,:) = [i, ii, B, SV, Z_num, Mr, Mo, FS_o, check_o, q_ult_uw, q_max, FS_bc, check_bc, Fh_r, sum_Fh, FS_s, check_s, check_r, check_p, overall_check, overallCost];
        indiv_geogrid_Dat{i,ii} = indiv_geogrid;
        rowN = rowN + 1;

        if overall_check == 1 
            workingsummaryData(rowNN,:) = [i, ii, B, SV, Z_num, Mr, Mo, FS_o, check_o, q_ult_uw, q_max, FS_bc, check_bc, Fh_r, sum_Fh, FS_s, check_s, check_r, check_p, overall_check, overallCost];
            rowNN = rowNN + 1;
        end
    end
end
sortedworkingsummaryData = sortrows(workingsummaryData,21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlmwrite('summary data of all results.csv', allsummaryData, 'delimiter', ',','precision',significantFigure);
dlmwrite('summary data of all working results_sorted.csv', sortedworkingsummaryData, 'delimiter', ',','precision',significantFigure);
dlmwrite('summary data of internal stability analysis - critical.csv', indiv_geogrid_Dat{sortedworkingsummaryData(1,1),sortedworkingsummaryData(1,2)}, 'delimiter', ',','precision',significantFigure);
disp('done')