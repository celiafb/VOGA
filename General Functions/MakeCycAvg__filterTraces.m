%% Filter position and velocity
function [filt,Data_calc,LE_V,RE_V,Data_cal,Data_In,LE_P,RE_P] = MakeCycAvg__filterTraces(type,filt_params_p,filt_params_v,te,ts,Data,keep_inds)
    %Make sure to add additional filter types to the filter parameter check
    %function
    %Position Filter Parameters
    filt.pos.l_x.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_p(1),'med');
    filt.pos.r_x.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_p(2),'med');
    filt.pos.l_y.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_p(3),'med');
    filt.pos.r_y.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_p(4),'med');
    filt.pos.l_z.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_p(5),'med');
    filt.pos.r_z.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_p(6),'med');
    filt.pos.l_x.spline = MakeCycAvg__filterParameterCheck(filt_params_p(7),'spline');
    filt.pos.r_x.spline = MakeCycAvg__filterParameterCheck(filt_params_p(8),'spline');
    filt.pos.l_y.spline = MakeCycAvg__filterParameterCheck(filt_params_p(9),'spline');
    filt.pos.r_y.spline = MakeCycAvg__filterParameterCheck(filt_params_p(10),'spline');
    filt.pos.l_z.spline = MakeCycAvg__filterParameterCheck(filt_params_p(11),'spline');
    filt.pos.r_z.spline = MakeCycAvg__filterParameterCheck(filt_params_p(12),'spline');
    filt.pos.l_x.sgolay = MakeCycAvg__filterParameterCheck(filt_params_p([13,19]),'sgolay');
    filt.pos.r_x.sgolay = MakeCycAvg__filterParameterCheck(filt_params_p([14,20]),'sgolay');
    filt.pos.l_y.sgolay = MakeCycAvg__filterParameterCheck(filt_params_p([15,21]),'sgolay');
    filt.pos.r_y.sgolay = MakeCycAvg__filterParameterCheck(filt_params_p([16,22]),'sgolay');
    filt.pos.l_z.sgolay = MakeCycAvg__filterParameterCheck(filt_params_p([17,23]),'sgolay');
    filt.pos.r_z.sgolay = MakeCycAvg__filterParameterCheck(filt_params_p([18,24]),'sgolay');
    %Velocity Filter Parameters
    filt.vel.irlssmooth = MakeCycAvg__filterParameterCheck(filt_params_v(1),'smooth');
    filt.vel.spline = MakeCycAvg__filterParameterCheck(filt_params_v(2),'spline');
    filt.vel.accelcutoff = MakeCycAvg__filterParameterCheck(filt_params_v(3),'acc');
    filt.vel.medianfilt = MakeCycAvg__filterParameterCheck(filt_params_v(4),'med');
    %Filter position
    LE_X = sgolay_filt(med_filt(Data.LE_Position_X,filt.pos.l_x.medianfilt),filt.pos.l_x.sgolay);
    RE_X = sgolay_filt(med_filt(Data.RE_Position_X,filt.pos.r_x.medianfilt),filt.pos.r_x.sgolay);
    LE_Y = sgolay_filt(med_filt(Data.LE_Position_Y,filt.pos.l_y.medianfilt),filt.pos.l_y.sgolay);
    RE_Y = sgolay_filt(med_filt(Data.RE_Position_Y,filt.pos.r_y.medianfilt),filt.pos.r_y.sgolay);
    LE_Z = sgolay_filt(med_filt(Data.LE_Position_Z,filt.pos.l_z.medianfilt),filt.pos.l_z.sgolay);
    RE_Z = sgolay_filt(med_filt(Data.RE_Position_Z,filt.pos.r_z.medianfilt),filt.pos.r_z.sgolay);
    % Recalculate Angular Velocity
    Data_In.Data_LE_Pos_X = spline_filt(te,LE_X,ts,filt.pos.l_x.spline);
    Data_In.Data_LE_Pos_Y = spline_filt(te,LE_Y,ts,filt.pos.l_y.spline);
    Data_In.Data_LE_Pos_Z = spline_filt(te,LE_Z,ts,filt.pos.l_z.spline);
    Data_In.Data_RE_Pos_X = spline_filt(te,RE_X,ts,filt.pos.r_x.spline);
    Data_In.Data_RE_Pos_Y = spline_filt(te,RE_Y,ts,filt.pos.r_y.spline);
    Data_In.Data_RE_Pos_Z = spline_filt(te,RE_Z,ts,filt.pos.r_z.spline);
    Data_In.Fs = Data.Fs;
    %First param says no initial coordinate transforms, second is a 
    %struct with the filtered position traces.
    Data_cal = angpos2angvel(1,Data_In);  
    %Filter velocity
    switch type
        case 1
            Data_calc.LE_Vel_X = spline_filt(ts,irls_filt(Data_cal.LE_Vel_X,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_X = spline_filt(ts,irls_filt(Data_cal.RE_Vel_X,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_Y = spline_filt(ts,irls_filt(Data_cal.LE_Vel_Y,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_Y = spline_filt(ts,irls_filt(Data_cal.RE_Vel_Y,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_Z = spline_filt(ts,irls_filt(Data_cal.LE_Vel_Z,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_Z = spline_filt(ts,irls_filt(Data_cal.RE_Vel_Z,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_LARP = spline_filt(ts,irls_filt(Data_cal.LE_Vel_LARP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_LARP = spline_filt(ts,irls_filt(Data_cal.RE_Vel_LARP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.LE_Vel_RALP = spline_filt(ts,irls_filt(Data_cal.LE_Vel_RALP,filt.vel.irlssmooth),ts,filt.vel.spline);
            Data_calc.RE_Vel_RALP = spline_filt(ts,irls_filt(Data_cal.RE_Vel_RALP,filt.vel.irlssmooth),ts,filt.vel.spline); 
        case 2
            %Turn points into dots and desaccade that way
            Data_calc = Data_cal;
            if ~isempty(filt.vel.accelcutoff)
                accelcutoff = filt.vel.accelcutoff;
                %Left Eye QPR
                d_LZ = [0;diff(Data_calc.LE_Vel_Z)];
                L_up_sacc = find(d_LZ>accelcutoff);
                L_up_sacc_s = L_up_sacc([true;diff(L_up_sacc) > 1]);
                for i = 1:length(L_up_sacc_s)
                    Data_calc.LE_Vel_Z(find(d_LZ(1:L_up_sacc_s(i)-1)<0,1,'last'):find(d_LZ(L_up_sacc_s(i)+1:end)<0,1,'first')+L_up_sacc_s(i)) = NaN;
                end
                L_down_sacc = find(d_LZ<-accelcutoff);
                L_down_sacc_s = L_down_sacc([true;diff(L_down_sacc) > 1]);
                for i = 1:length(L_down_sacc_s)
                    Data_calc.LE_Vel_Z(find(d_LZ(1:L_down_sacc_s(i)-1)>0,1,'last'):find(d_LZ(L_down_sacc_s(i)+1:end)>0,1,'first')+L_down_sacc_s(i)) = NaN;
                end
                LZ_nan = find(isnan(Data_calc.LE_Vel_Z));
                small_spac = find(diff(LZ_nan) > 1 & diff(LZ_nan) < 3);
                for i = 1:length(small_spac)
                    Data_calc.LE_Vel_Z(LZ_nan(small_spac(i)):LZ_nan(small_spac(i)+1)) = NaN;
                end
                Data_calc.LE_Vel_LARP(isnan(Data_calc.LE_Vel_Z)) = NaN;
                Data_calc.LE_Vel_RALP(isnan(Data_calc.LE_Vel_Z)) = NaN;
                Data_calc.LE_Vel_X(isnan(Data_calc.LE_Vel_Z)) = NaN;
                Data_calc.LE_Vel_Y(isnan(Data_calc.LE_Vel_Z)) = NaN;
                %Right eye QPR
                d_RZ = [0;diff(Data_calc.RE_Vel_Z)];
                R_up_sacc = find(d_RZ>accelcutoff);
                R_up_sacc_s = R_up_sacc([true;diff(R_up_sacc) > 1]);
                for i = 1:length(R_up_sacc_s)
                    Data_calc.RE_Vel_Z(find(d_RZ(1:R_up_sacc_s(i)-1)<0,1,'last'):find(d_RZ(R_up_sacc_s(i)+1:end)<0,1,'first')+R_up_sacc_s(i)) = NaN;
                end
                R_down_sacc = find(d_RZ<-accelcutoff);
                R_down_sacc_s = R_down_sacc([true;diff(R_down_sacc) > 1]);
                for i = 1:length(R_down_sacc_s)
                    Data_calc.RE_Vel_Z(find(d_RZ(1:R_down_sacc_s(i)-1)>0,1,'last'):find(d_RZ(R_down_sacc_s(i)+1:end)>0,1,'first')+R_down_sacc_s(i)) = NaN;
                end
                RZ_nan = find(isnan(Data_calc.RE_Vel_Z));
                small_spac = find(diff(RZ_nan) > 1 & diff(RZ_nan) < 3);
                for i = 1:length(small_spac)
                    Data_calc.LE_Vel_Z(RZ_nan(small_spac(i)):RZ_nan(small_spac(i)+1)) = NaN;
                end
                Data_calc.RE_Vel_LARP(isnan(Data_calc.RE_Vel_Z)) = NaN;
                Data_calc.RE_Vel_RALP(isnan(Data_calc.RE_Vel_Z)) = NaN;
                Data_calc.RE_Vel_X(isnan(Data_calc.RE_Vel_Z)) = NaN;
                Data_calc.RE_Vel_Y(isnan(Data_calc.RE_Vel_Z)) = NaN; 
            end                        
            Data_calc.LE_Vel_X = med_filt(Data_calc.LE_Vel_X,filt.vel.medianfilt);
            Data_calc.LE_Vel_Y = med_filt(Data_calc.LE_Vel_Y,filt.vel.medianfilt);
            Data_calc.LE_Vel_Z = med_filt(Data_calc.LE_Vel_Z,filt.vel.medianfilt);
            Data_calc.LE_Vel_LARP = med_filt(Data_calc.LE_Vel_LARP,filt.vel.medianfilt);
            Data_calc.LE_Vel_RALP = med_filt(Data_calc.LE_Vel_RALP,filt.vel.medianfilt);
            Data_calc.RE_Vel_X = med_filt(Data_calc.RE_Vel_X,filt.vel.medianfilt);
            Data_calc.RE_Vel_Y = med_filt(Data_calc.RE_Vel_Y,filt.vel.medianfilt);
            Data_calc.RE_Vel_Z = med_filt(Data_calc.RE_Vel_Z,filt.vel.medianfilt);
            Data_calc.RE_Vel_LARP = med_filt(Data_calc.RE_Vel_LARP,filt.vel.medianfilt);
            Data_calc.RE_Vel_RALP = med_filt(Data_calc.RE_Vel_RALP,filt.vel.medianfilt);
    end
    %Initialize vectors with the cycles
    LE_V.LHRH = Data_calc.LE_Vel_Z(keep_inds);
    LE_V.LARP = Data_calc.LE_Vel_LARP(keep_inds);
    LE_V.RALP = Data_calc.LE_Vel_RALP(keep_inds);
    LE_V.X = Data_calc.LE_Vel_X(keep_inds);
    LE_V.Y = Data_calc.LE_Vel_Y(keep_inds);
    RE_V.LHRH = Data_calc.RE_Vel_Z(keep_inds);
    RE_V.LARP = Data_calc.RE_Vel_LARP(keep_inds);
    RE_V.RALP = Data_calc.RE_Vel_RALP(keep_inds);
    RE_V.X = Data_calc.RE_Vel_X(keep_inds);
    RE_V.Y = Data_calc.RE_Vel_Y(keep_inds);
    
    %And for Position MRC 3/2/21
    LE_P.LHRH = Data_In.Data_LE_Pos_Z(keep_inds);
    LE_P.X = Data_In.Data_LE_Pos_X(keep_inds);
    LE_P.Y = Data_In.Data_LE_Pos_Y(keep_inds);
    RE_P.LHRH = Data_In.Data_RE_Pos_Z(keep_inds);
    RE_P.X = Data_In.Data_RE_Pos_X(keep_inds);
    RE_P.Y = Data_In.Data_RE_Pos_Y(keep_inds);
end