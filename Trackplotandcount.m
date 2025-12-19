function [phasestandarddeviation,dopstandarddeviation,CN0MEAN]=Trackplotandcount(runtime,carrerr,DOPERRORre,CN0L,FCODE,FCARR,codeerr)
%% 本函数用于统计跟踪误差
tracktime=runtime*1000-10;
limtimems=0.5*1000;%收敛时间
%数据截断：
adjusted_phase = (carrerr - round(carrerr));

% adjusted_phase = (carrerr);
% for i=1:1:tracktime
%     if(adjusted_phase(i)>0.5)
%         adjusted_phase(i)=adjusted_phase(i)-1;
%     elseif(adjusted_phase(i)<-0.5)
%          adjusted_phase(i)=adjusted_phase(i)+1;
%    
%     end
%     
% end
for i=1:1:tracktime
    if(adjusted_phase(i)>0.25)
        adjusted_phase(i)=adjusted_phase(i)-0.5;
    elseif(adjusted_phase(i)<-0.25)
        adjusted_phase(i)=adjusted_phase(i)+0.5;
    end
end
adjusted_phase=adjusted_phase(limtimems:1:tracktime);
adjusted_phase_code=(codeerr - round(codeerr));
adjusted_phase_code=adjusted_phase_code(limtimems:1:tracktime);
FCODEresult=FCODE(limtimems:1:tracktime);
FCARRresult=FCARR(limtimems:1:tracktime);
DOPERRORreresult=DOPERRORre(limtimems:1:tracktime);
CN0Lresult=CN0L(limtimems:1:tracktime);

% figure,plot(FCODEresult),title('码多普勒'),xlabel("ms");
% figure,plot(FCARRresult),title('载波多普勒'),xlabel("2ms");
% figure,plot(DOPERRORreresult),title('本地与接受多普勒差值 / Hz'),xlabel("ms");
% figure,plot(adjusted_phase),title('载波相位误差 / 周期数'),xlabel("ms");
% figure,plot(adjusted_phase_code),title('码相位误差/ chip'),xlabel("ms");
% figure,plot(abs(CN0Lresult)),title('载噪比估计')

% 载波相位误差计算：
adjusted_phasemean_value=mean(adjusted_phase);
% 载波相位均方误差
mse = mean((adjusted_phase - adjusted_phasemean_value).^2);
% standard deviation
phasestandarddeviation=sqrt(mse);

%载波多普勒误差计算：
DOPERRORremean_value=mean(DOPERRORreresult);
%均方误差
dopmse = mean((DOPERRORremean_value - DOPERRORreresult).^2);
%standard deviation
dopstandarddeviation=sqrt(dopmse);


%载噪比估计
CN0MEAN=mean(abs(CN0Lresult));

% %LKF  PD的观测Rk方差  
% sigma_meas = std(PD);
% fprintf('实际观测 pd 的 RMS = %.4f 周期\n', sigma_meas);
% %LKF设置的R_KF
% R_KFRMS=sqrt(R_kf);
% %载波相位180°模糊时用下面的
% standarddeviation_adjusted_phase_half=sqrt(mean((carrerr - mean(carrerr)).^2));
fprintf('载波相位RMSE(cycles) = %.8f  多普勒RMSE(Hz) = %.8f   载噪比均值 = %.8f\n', phasestandarddeviation,dopstandarddeviation,CN0MEAN);
end
