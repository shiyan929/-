
%      [ 莫来石, t_p ] | [ 微间隙, d ] | [ 莫来石, t_p ]
%
%   1) 稀薄气体导热计算
%        lambda  = kB*Tf / ( sqrt(2)*pi*d_mol^2*P )     平均自由程
%        Kn      = lambda / d          克努森数
%        beta    = (2-alpha)/alpha * (2*gamma/(gamma+1)) / Pr   温度跳跃系数
%        k_g,eff = k_g,0 / (1 + 2*beta*Kn)    有效气体导热
%        h_gas   = k_g,eff / d  =  k_g,0 / (d + 2*beta*lambda)  气体间隙导热系数
%
%   2) 黑体辐射
%        eps_eff = 1/(1/eps_h + 1/eps_c - 1)
%        h_rad   = eps_eff*sigma*(Th+Tc)*(Th^2+Tc^2)
%
%   3) 固体桥柱
%        h_pil   = k_pillar*phi / d
%
%  总间隙：
%        h_gap   = h_gas + h_rad + h_pil   间隙导热系数
%        R_total = 2*t_p/k_m + 1/h_gap    总热阻
%        L_total = 2*t_p + d         总厚度
%        k_stack = L_total / R_total   等效导热系数

clear; close all; clc;

here   = pwd;
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% 常量定义
C.kB    = 1.380649e-23;    % 玻尔兹曼常数        [J/K]
C.sigma = 5.670374419e-8;  % 斯特藩-玻尔兹曼常数 [W/(m^2*K^4)]
C.p_atm = 101325;          % 标准大气压       [Pa]

%% 输入定义
D.Th        = 1073.15;   % 热端温度            [K]  (~800 C)
D.Tc        = 293.15;    % 冷端温度            [K]  (~20 C)
D.t_plate   = 2.0e-3;    % 莫来石板厚度        [m]
D.k_mullite = 4.0;       % 莫来石板导热系数    [W/(m*K)]
D.d_mol     = 3.70e-10;  % 有效分子直径，空气  [m]
D.gamma     = 1.4;       % 比热容，空气       
D.Pr        = 0.71;      % 空气普朗特数
D.alpha     = 0.90;      % 热适应系数
D.emis_h    = 0.80;      % 表面发射率，热端
D.emis_c    = 0.80;      % 表面发射率，冷端
D.k_pillar  = 4.0;       % 桥柱材料导热系数
D.phi       = 0.00;      % 桥柱面积
%% 辐射属性
Tf      = 0.5*(D.Th + D.Tc);
kg0     = air_k(Tf);     %计算空气导热系数
beta    = (2-D.alpha)/D.alpha * (2*D.gamma/(D.gamma+1)) / D.Pr;   % 温度跳跃系数
eps_eff = 1/(1/D.emis_h + 1/D.emis_c - 1);     % 两平行板的有效发射率
h_rad   = eps_eff * C.sigma * (D.Th+D.Tc)*(D.Th^2+D.Tc^2);   % 辐射导热系数
lam_atm = C.kB*Tf/(sqrt(2)*pi*D.d_mol^2*C.p_atm);   %标准大气压下的平均自由程

fprintf('膜温度             Tf      = %9.2f K\n',  Tf);
fprintf('空气本体导热系数   kg0     = %9.5f W/(m*K)\n', kg0);
fprintf('温度跳跃系数       beta    = %9.4f\n', beta);
fprintf('平均自由程 @1atm   lam     = %9.2f nm\n', lam_atm*1e9);
fprintf('有效发射率         eps_eff = %9.4f\n', eps_eff);
fprintf('辐射导热系数       h_rad   = %9.3f W/(m^2*K)\n', h_rad);
fprintf('总热阻（双板）     R_p     = %9.3e m^2*K/W\n', ...
        2*D.t_plate/D.k_mullite);
fprintf('莫来石导热系数     k       = %9.3f W/(m*K)\n\n', ...
        D.k_mullite);

%% ========================================================================
%  FIG 1 稀薄气体通用曲线 k_gas_eff/k_gas0 vs Kn
%% ========================================================================
Kn_v  = logspace(-3, 2, 400);   %从0.001~100，生成400个对数等分点
ratio = 1 ./ (1 + 2*beta*Kn_v);     %计算有效气体导热

f1 = figure('Visible','off','Position',[80 80 820 560]);
semilogx(Kn_v, ratio, 'LineWidth', 2); grid on; hold on;
yline(0.5, ':', 'Color', [.5 .5 .5]);
xline(0.01,'--','Color',[.45 .45 .45]);
xline(0.1 ,'--','Color',[.45 .45 .45]);
xline(10  ,'--','Color',[.45 .45 .45]);
xlim([1e-3 1e2]); ylim([0 1.05]);
xlabel('Knudsen number   Kn = \lambda / d');
ylabel('k_{g,eff} / k_{g,0}');
title(sprintf('Rarefaction suppression of gap conductivity (\\alpha = %.2f)', D.alpha));
text(0.02, 0.06, 'slip', 'FontSize', 9);
text(0.25, 0.06, 'transition', 'FontSize', 9);
text(2.0 , 0.06, 'free molecular', 'FontSize', 9);
exportgraphics(f1, fullfile(outdir,'fig1_rarefaction_curve.png'), 'Resolution',150);
close(f1);

%% ========================================================================
%  FIG 2 -- k_stack（等效导热系数） vs 不同压力下，不同间隙厚度
%% ========================================================================
d_vec = logspace(-7, -2, 400);                 % 0.1 um .. 10 mm
P_set = [100 1000 1e4 C.p_atm];                % Pa
P_lab = {'100 Pa','1000 Pa','10 kPa','1 atm'};
col   = lines(numel(P_set));

f2 = figure('Visible','off','Position',[80 80 900 600]);
hold on; grid on; box on;
for i = 1:numel(P_set)
    r = model_stack(d_vec, P_set(i), D, C, kg0, beta, h_rad);
    semilogx(d_vec*1e6, r.k_stack, 'LineWidth', 2, 'Color', col(i,:), ...
             'DisplayName', P_lab{i});
end
yline(D.k_mullite, 'k--', 'DisplayName', 'bare mullite');
set(gca,'YScale','log');
xlabel('gap thickness   d  [\mum]');
ylabel('stack effective conductivity   k_{stack}  [W/(m*K)]');
title('Effective conductivity of the plate/gap/plate stack');
legend('Location','northeast'); xlim([0.1 1e4]);
exportgraphics(f2, fullfile(outdir,'fig2_kstack_vs_gap.png'), 'Resolution',150);
close(f2);

%% ========================================================================
%  FIG 3 -- 等效热导系数和Kn等值线
%% ========================================================================
d2 = logspace(-7, -2, 240);   %间隙厚度步长
P2 = logspace(0, 5.2, 200);   %压力步长
[DD, PP] = meshgrid(d2, P2);
R2 = model_stack(DD, PP, D, C, kg0, beta, h_rad);    %对每个点进行计算

f3 = figure('Visible','off','Position',[80 60 980 720]);
contourf(log10(DD), log10(PP), log10(R2.k_stack), 28, 'LineColor','none');
hold on;
[c,h] = contour(log10(DD), log10(PP), R2.Kn, log10([0.1 1 10]), ...
                'k', 'LineWidth', 1.6, 'ShowText','on');
clabel(c, h, 'FontSize', 9, 'Color','k');
xt = -7:1:-2;
set(gca,'XTick',xt,'XTickLabel', ...
    arrayfun(@(v) sprintf('10^{%d}',v), xt, 'UniformOutput', false));
yt = 0:1:5;
set(gca,'YTick',yt,'YTickLabel', ...
    arrayfun(@(v) sprintf('10^{%d}',v), yt, 'UniformOutput', false));
xlabel('gap thickness   d  [m]');
ylabel('gas pressure   P  [Pa]');
title('Design map: log_{10} k_{stack}  with Kn contours');
cb = colorbar; cb.Label.String = 'log_{10}( k_{stack} )  [W/(m*K)]';
exportgraphics(f3, fullfile(outdir,'fig3_design_map.png'), 'Resolution',150);
close(f3);

%% ========================================================================
%  FIG 4 -- 传热机制分析
%% ========================================================================
P_fix = 1000;    %恒定压力
rb = model_stack(d_vec, P_fix, D, C, kg0, beta, h_rad);

f4 = figure('Visible','off','Position',[80 80 900 600]);
loglog(d_vec*1e6, rb.h_gas, 'LineWidth', 2); hold on; grid on;
loglog(d_vec*1e6, rb.h_rad, 'LineWidth', 2);
loglog(d_vec*1e6, rb.h_gap, 'k--', 'LineWidth', 2);
set(gca,'YScale','log');
xlabel('gap thickness   d  [\mum]');
ylabel('conductance   h  [W/(m^2*K)]');
title(sprintf('Heat-transfer mechanism breakdown at P = %d Pa', P_fix));
legend('gas (rarefied)','radiation','total gap', 'Location','best');
exportgraphics(f4, fullfile(outdir,'fig4_mechanism_breakdown.png'), 'Resolution',150);
close(f4);

%% ========================================================================
%  FIG 5 -- 莫来石厚度敏感性
%% ========================================================================
t_set = logspace(-4, -2.3, 40);               %对数点设计
d_des = 20e-6;                                 %设计间隙
P_des = 1000;                                  % 压力设计
k_t = zeros(size(t_set));
for i = 1:numel(t_set)
    Dz = D; Dz.t_plate = t_set(i);
    r  = model_stack(d_des, P_des, Dz, C, kg0, beta, h_rad);
    k_t(i) = r.k_stack;
end

f5 = figure('Visible','off','Position',[80 80 860 570]);
semilogx(t_set*1e3, k_t, 'LineWidth', 2); grid on;
xlabel('single plate thickness   t_p  [mm]');
ylabel('stack effective conductivity   k_{stack}  [W/(m*K)]');
title(sprintf('Plate thickness sensitivity (d = %g \\mum, P = %d Pa)', ...
      d_des*1e6, P_des));
exportgraphics(f5, fullfile(outdir,'fig5_plate_thickness.png'), 'Resolution',150);
close(f5);

%% ========================================================================
%  DESIGN TABLE -- gap that puts Kn in the transition band
%% ========================================================================
P_tab = [1 10 100 1000 1e4 C.p_atm];
fprintf('========== DESIGN TABLE (transition regime) ==========\n');
fprintf('%10s %12s %12s %12s %12s\n', ...
        'P [Pa]', 'lam [m]', 'd(Kn=1)[m]', 'd(Kn=0.1)[m]', 'k_stack');
fprintf('%10s %12s %12s %12s %12s\n', ...
        '', '', '', '', '[W/mK]');

T = zeros(numel(P_tab), 4);
for i = 1:numel(P_tab)
    P = P_tab(i);
    lam = C.kB*Tf/(sqrt(2)*pi*D.d_mol^2*P);
    d1  = lam/1.0;
    d01 = lam/0.1;
    r   = model_stack(d1, P, D, C, kg0, beta, h_rad);
    T(i,:) = [P, lam, d1, r.k_stack];
    fprintf('%10.0f %12.3e %12.3e %12.3e %12.4f\n', P, lam, d1, d01, r.k_stack);
end
fprintf('\n');

writematrix(T, fullfile(outdir,'design_table.csv'));

%% export the 2D map for further analysis
d2o = repmat(d2, numel(P2), 1);
P2o = repmat(P2(:), 1, numel(d2));
mapOut = [d2o(:), P2o(:), R2.Kn(:), R2.k_stack(:)];
writematrix(mapOut, fullfile(outdir,'design_map.csv'));

fprintf('Figures + CSV written to: %s\n', outdir);
fprintf('================ DONE ================\n');

%% ========================================================================
%  保存图片
%% ========================================================================
function k = air_k(T)
%AIR_K  Bulk thermal conductivity of air, power-law fit  [W/(m*K)]
k = 0.02624 * (T/300).^0.8646;
end

function r = model_stack(d, P, D, C, kg0, beta, h_rad)
%MODEL_STACK  Evaluate the plate/gap/plate stack.
%   d [m] gap thickness (scalar or array), P [Pa] pressure.

Tf    = 0.5*(D.Th + D.Tc);
lam   = C.kB*Tf ./ (sqrt(2)*pi*D.d_mol^2 .* P);
Kn    = lam ./ d;
kg    = kg0 ./ (1 + 2*beta.*Kn);
h_gas = kg ./ d;
h_pil = D.k_pillar*D.phi ./ d;
h_gap = h_gas + h_rad + h_pil;
Rgap  = 1 ./ h_gap;
Rtot  = 2*D.t_plate/D.k_mullite + Rgap;
Ltot  = 2*D.t_plate + d;
k_st  = Ltot ./ Rtot;

r = struct( ...
    'lam',     lam, ...
    'Kn',      Kn, ...
    'kg',      kg, ...
    'h_gas',   h_gas, ...
    'h_pil',   h_pil, ...
    'h_rad',   h_rad*ones(size(d)), ...
    'h_gap',   h_gap, ...
    'k_stack', k_st, ...
    'k_layer', h_gap .* d, ...
    'Rtot',    Rtot, ...
    'Ltot',    Ltot);
end
