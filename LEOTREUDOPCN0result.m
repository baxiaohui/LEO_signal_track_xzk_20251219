clear; close all; clc;

%% ================= 原始数据（已清理可直接运行） =================
CHR = {
[45 44 43 42 41 40 39 38 37 36 ];
[45 44 43 42 41 40 39 38 37 36 35 ];
[45 44 43 42 41 40 39 38 37 36 35 34 33 ];
[45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 ];
[45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26]
};

CarrierPhaseRMSE = {
[    0.008064283 0.008957145 0.009975292 0.011143741 0.012470525 0.013997535 0.015754716 0.017892865 0.020772097 0.027228778 ];
[0.006499444	0.007312224	0.008231813	0.009285317	0.010483572	0.011885213	0.013521811	0.015440533	0.01778347	0.020661687	0.024073585
];
[0.004840653	0.005446675	0.006132969	0.006914708	0.007802529	0.010231497	0.011644469	0.01328732	0.015223438	0.015260508	0.017725988	0.024017904	0.029243714
];
[0.005735104	0.006350863	0.007026373	0.007779829	0.008605906	0.009544479	0.010580326	0.011720887	0.013008413	0.014469058	0.016029806	0.017786401	0.019874772	0.022207091	0.024871415	0.028175411	0.032117603	0.037443155
];
[0.004903301	0.005420753	0.005995292	0.006633092	0.007330751	0.008117441	0.008991984	0.009967888	0.011040673	0.01225381	0.013589028	0.015052785	0.016781878	0.018670995	0.020861419	0.023465351	0.026417054	0.029971666	0.034009125	0.039493473
];

};

DopplerRMSE = {
[1.151741528	1.27828062	1.42276309	1.588049397	1.776366848	1.992616493	2.240281283	2.535891713	2.918732001	3.57 ];
[0.24136907	0.259722925	0.281302859	0.307006055	0.337190696	0.373344178	0.416307704	0.467614954	0.532198114	0.611786182	0.714955235
];
[0.162130452	0.16422454	0.166852844	0.170170987	0.174333655	0.210588007	0.224537944	0.241798306	0.263320912	0.272027848	0.28271639	0.370215809	0.559039063
];
[0.200275969	0.206683274	0.213696487	0.221640156	0.230391156	0.240443066	0.251456278	0.26354029	0.277381333	0.293101915	0.309450707	0.3277185	0.349781375	0.374322253	0.401560005	0.436807169	0.478827489	0.538637105
];
[0.163184754	0.164542119	0.166099182	0.167886875	0.169901987	0.17226258	0.174968321	0.178088737	0.181629975	0.185774745	0.190408387	0.195622661	0.201973114	0.208977143	0.217329666	0.227683257	0.239519962	0.254229452	0.270846124	0.294959363
];

};

schemeNames = {'PLL','LKF','DS-LKF','EKF','DS-EKF'};

%% ================= 样式设置（色盲安全 + 期刊风格） =================
pal = [ ...                     % Okabe–Ito palette
0.00 0.45 0.70;   % 蓝   PLL
0.90 0.60 0.00;   % 橙   LKF
0.00 0.62 0.45;   % 青   DS-LKF
0.80 0.47 0.65;   % 品红 EKF
0.35 0.70 0.90];  % 天蓝 DS-EKF
markers    = {'o','s','d','^','v'};
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
% plot_ber('Coherent BER',      BERCoherent,     CHR, schemeNames, pal, markers, linestyles, lw, ms, mkstep, 1e-6,  'best');
% plot_ber('Noncoherent BER',   BERNonCoherent,  CHR, schemeNames, pal, markers, linestyles, lw, ms, mkstep, 1e-6,  'southeast');
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
xlim([26 45]); xticks(26:45); ylim([0 0.10]);
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
xlim([26 45]); xticks(26:45); ylim([0 4]);
legend('Location','northeast','NumColumns',3,'Box','off');
title('(b) Doppler RMSE');
set(findall(fig,'Type','axes'),'LooseInset',[0 0 0 0]);
exportgraphics(fig,'LEOTREUDOPCN0result.png','Resolution',600);
exportgraphics(fig,'LEOTREUDOPCN0result.pdf','ContentType','vector');
exportgraphics(fig,'LEOTREUDOPCN0result.emf','ContentType','vector');
% ====== 高质量矢量导出 ======
% exportgraphics(fig,'RMSE_vs_CN0.pdf','ContentType','vector');      % PDF（推荐）
% print(fig,'RMSE_vs_CN0.eps','-depsc','-painters');               % EPS（投稿常用）
% print(fig,'RMSE_vs_CN0.emf','-dmeta');                           % EMF（Word/PPT）
