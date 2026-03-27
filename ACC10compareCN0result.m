clear; close all; clc;

%% ================= 原始数据（已清理可直接运行） =================
CHR = {
[45 44 43 42 41 40 39 38 37 36 35];
[45 44 43 42 41 40 39 38 37 36 35 34 33];
[45 44 43 42 41 40 39 38 37 36 35 34 33 32];
[45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24];
[45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23]
};

BERCoherent = [ % 每行 11 个值（对应 CHR{i} 的前 11 个点）
NaN 1.04E-05 1.66E-04 1.20E-03 4.93E-03 1.60E-02 4.00E-02 7.62E-02 1.65E-01 2.27E-01 2.87E-01;
NaN 1.04E-05 1.76E-04 1.25E-03 4.84E-03 1.61E-02 3.97E-02 7.62E-02 1.25E-01 1.80E-01 2.36E-01;
NaN 1.04E-05 1.76E-04 1.22E-03 4.82E-03 1.60E-02 3.96E-02 7.59E-02 1.25E-01 1.79E-01 2.35E-01;
NaN 1.04E-05 1.87E-04 1.24E-03 5.44E-03 1.60E-02 3.96E-02 7.76E-02 1.27E-01 1.79E-01 2.34E-01;
NaN 1.04E-05 1.76E-04 1.19E-03 4.81E-03 1.60E-02 3.95E-02 7.58E-02 1.24E-01 1.78E-01 2.34E-01
];

BERNonCoherent = [
NaN 8.29E-05 8.19E-04 4.42E-03 1.67E-02 4.40E-02 8.96E-02 1.49E-01 2.14E-01 2.75E-01 3.29E-01;
NaN 7.26E-05 8.19E-04 4.48E-03 1.66E-02 4.42E-02 8.96E-02 1.49E-01 2.14E-01 2.76E-01 3.29E-01;
NaN 8.29E-05 8.19E-04 4.47E-03 1.66E-02 4.40E-02 8.97E-02 1.49E-01 2.14E-01 2.75E-01 3.29E-01;
1.04E-05 8.29E-05 7.98E-04 4.19E-03 1.68E-02 4.40E-02 8.97E-02 1.50E-01 2.16E-01 2.77E-01 3.29E-01;
NaN 7.26E-05 7.88E-04 4.48E-03 1.66E-02 4.41E-02 8.97E-02 1.49E-01 2.14E-01 2.76E-01 3.29E-01
];

BERLDPCCoherent = [
NaN NaN NaN NaN NaN NaN NaN NaN 3.92E-03 1.01E-01 2.36E-01;
NaN NaN NaN NaN NaN NaN NaN NaN 5.18E-05 4.98E-03 1.04E-01;
NaN NaN NaN NaN NaN NaN NaN NaN 4.15E-05 4.63E-03 1.02E-01;
NaN NaN NaN NaN NaN NaN NaN NaN 2.07E-05 5.83E-03 9.93E-02;
NaN NaN NaN NaN NaN NaN NaN NaN 2.07E-05 4.42E-03 9.74E-02
];

CarrierPhaseRMSE = {
[0.006515601 0.007303839 0.008205618 0.009223823 0.010392183 0.011707163 0.013208556 0.014899917 0.017924979 0.021687776 0.028841124];
[0.005572679 0.00625071 0.00703577 0.007920687 0.008956744 0.0101369 0.011525718 0.013163769 0.015063923 0.017227635 0.020535813 0.024028712 0.029378407];
[0.004992135 0.00560272 0.006310523 0.00710573 0.008029706 0.009087994 0.010321836 0.011816943 0.013525315 0.015445002 0.018065521 0.020955733 0.024472042 0.028891488];
[0.005258401 0.005798043 0.006404491 0.007057937 0.007787085 0.008620647 0.009531863 0.01053925 0.011665127 0.012926624 0.014299453 0.015771 0.0175548 0.019395005 0.021623265 0.023991005 0.026803889 0.030219068 0.034403456 0.040100002 0.046730555 0.057923455];
[0.00459423 0.005069746 0.005607063 0.006193673 0.00684202 0.007587145 0.008384791 0.00929365 0.010310253 0.011415433 0.012628883 0.013943272 0.015508479 0.017247744 0.019107468 0.021165902 0.023619411 0.02630674 0.029440529 0.032821982 0.036796727 0.040346526 0.047616692]
};

DopplerRMSE = {
[0.832012688 0.932286513 1.045321322 1.174487641 1.320649988 1.487078563 1.676342278 1.892002267 2.212064062 2.603127972 3.281391524];
[0.1255095 0.139204959 0.155271357 0.173458015 0.195092667 0.219622369 0.248862254 0.28344015 0.323503304 0.368449098 0.440133024 0.5153668 0.633599435];
[0.070520975 0.076358441 0.083364202 0.091440786 0.101074243 0.112291545 0.125680096 0.141331028 0.160269811 0.181238376 0.210633824 0.243740705 0.284104711 0.337895454];
[0.1043236 0.110608807 0.11766307 0.124978348 0.132998359 0.142120465 0.151884708 0.162336338 0.173919642 0.186603173 0.199848535 0.213534981 0.230416167 0.246408911 0.266299843 0.286359351 0.31053045 0.340010838 0.375748258 0.427967649 0.485382534 0.779235208];
[0.060762271 0.063261675 0.066087138 0.069166091 0.072555575 0.076476289 0.080708775 0.085475449 0.090881758 0.096662773 0.102920933 0.109776969 0.1179127 0.126556094 0.135809647 0.145853466 0.158058711 0.171239884 0.186459245 0.203327617 0.221852026 0.234501984 0.269015038]
};

schemeNames = {'PLL','LKF','DS-LKF','EKF','DS-EKF'};

%% ================= 样式设置（色盲安全 + 期刊风格） =================
pal = [ ... % Okabe–Ito palette
0.00 0.45 0.70; % 蓝 PLL
0.90 0.60 0.00; % 橙 LKF
0.00 0.62 0.45; % 青 DS-LKF
0.80 0.47 0.65; % 品红 EKF
0.35 0.70 0.90]; % 天蓝 DS-EKF
markers = {'o','s','d','^','v'};
linestyles = {'-','--','-.',':','-'};
lw = 2; ms = 6;
mkstep = 1; % 数据点不多，全部显示

% 全局字体（Times + 粗体）
set(groot,'defaultAxesFontName','Times New Roman', ...
'defaultTextFontName','Times New Roman', ...
'defaultAxesFontSize',24, ...
'defaultTextFontSize',24, ...
'defaultAxesFontWeight','bold', ...
'defaultTextFontWeight','bold');

% %% ================= 绘制 BER（三张图） =================
% plot_ber('Coherent BER', BERCoherent, CHR, schemeNames, pal, markers, linestyles, lw, ms, mkstep, 1e-6, 'best');
% plot_ber('Noncoherent BER', BERNonCoherent, CHR, schemeNames, pal, markers, linestyles, lw, ms, mkstep, 1e-6, 'southeast');
% plot_ber('LDPC-Coherent BER', BERLDPCCoherent, CHR, schemeNames, pal, markers, linestyles, lw, ms, mkstep, 1e-10, 'southeast');

%% ================= 绘制 RMSE（上下 2×1 紧凑） =================
fig = figure('Color','w','Units','normalized','Position',[0.06 0.06 0.62 0.86]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% (a) Carrier Phase RMSE
nexttile; hold on; grid on;
for i=1:5
x = CHR{i}; y = CarrierPhaseRMSE{i};
plot(x, y, 'LineStyle',linestyles{i}, 'Color',pal(i,:), ...
'Marker',markers{i}, 'MarkerSize',ms, 'LineWidth',lw, ...
'MarkerIndices',1:mkstep:numel(x), 'DisplayName',schemeNames{i});
end
xlabel('C/N_0 (dBHz)'); ylabel('Carrier Phase RMSE (cycles)');
xlim([23 45]); xticks(23:45); ylim([0 0.10]);
legend('Location','northeast','NumColumns',3,'Box','off');
title('(a) Carrier Phase RMSE');

% (b) Doppler RMSE
nexttile; hold on; grid on;
for i=1:5
x = CHR{i}; y = DopplerRMSE{i};
plot(x, y, 'LineStyle',linestyles{i}, 'Color',pal(i,:), ...
'Marker',markers{i}, 'MarkerSize',ms, 'LineWidth',lw, ...
'MarkerIndices',1:mkstep:numel(x), 'DisplayName',schemeNames{i});
end
xlabel('C/N_0 (dBHz)'); ylabel('Doppler RMSE (Hz)');
xlim([23 45]); xticks(23:45); ylim([0 4]);
legend('Location','northeast','NumColumns',3,'Box','off');
title('(b) Doppler RMSE');



% 进一步去白边（导出时更明显）
set(findall(fig,'Type','axes'),'LooseInset',[0 0 0 0]);
% try
%     fig.WindowState = 'maximized';           % R2021a 支持
% catch
%     % 万一不支持，就用下面这一句
% end
% set(fig,'Units','normalized','OuterPosition',[0 0 1 1]); % 更通用、更稳
% 
% drawnow;  % 确保全屏生效后再导出

% ===== 导出（自动裁剪边缘白边）=====
exportgraphics(fig,'ACC10compareCN0result.png','Resolution',600);
exportgraphics(fig,'ACC10compareCN0result.pdf','ContentType','vector');
exportgraphics(fig,'ACC10compareCN0result.emf','ContentType','vector');
