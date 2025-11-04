% --- USTAWIENIA FIGURY ---
fig = figure('Color', 'w', ...
             'Units', 'pixels', ...
             'Position', [100, 100, 1300, 420]);
set(groot, 'defaultAxesFontSize', 20)
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

left_margins = [0.055, 0.38, 0.705];
width = 0.26;
height = 1;

% --- SUBPLOT 1 ---
ax1 = axes('Position', [left_margins(1), 0.055, width, height]);
eco_stacj_cf_WS_step([0.1 0.3 0.5], 4, 9999, 100);

% --- SUBPLOT 2 ---
ax2 = axes('Position', [left_margins(2), 0.055, width, height]);
eco_stacj_cf_WS_step([0.1 0.3 0.5], 4, 8, 100);

% --- SUBPLOT 3 ---
ax3 = axes('Position', [left_margins(3), 0.055, width, height]);
eco_stacj_cf_WS_step([0.1 0.3 0.5], 4, 8, 5);

% --- ZAPIS ---
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardcopy', 'off');
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches', 'PaperSize', [pos(3) pos(4)], ...
    'PaperPosition', [0 0 pos(3) pos(4)]);
% print(fig, 'eco_stacj_cf_WS_step_subplot.pdf', '-dpdf', '-r300')

% --- FUNKCJA ---
function eco_stacj_cf_WS_step(p, q, k, beta)
    % funkcja dla ustalonego p, q i k rysuje:
    % stacj. c(f) z MFA, PA i MC dla grafu WS
    % p mo¿e byæ wektorem

    points = ['^', 's', 'o', 'd', '*', 'h', 'v', 'x'];
    colors = [224, 159, 62; 158, 42, 43; 87, 117, 144]./255;  % yellow, red, blue
    f = 0:0.01:1;
    LH = [];
    leg = [];

    for pi = 1:length(p)
        % --- MFA ---
        c_st_mfa = [];
        dt = 0.25;
        T = dt:dt:800;
        for j = 1:length(f)
            c_mfa = 0;
            for t = 1:length(T)
                c_mfa = [c_mfa, c_mfa(end) + dt*(1 - p(pi))*((1 - c_mfa(end))*c_mfa(end)^q - c_mfa(end)*(1 - c_mfa(end))^q) ...
                    + dt*p(pi)*(f(j) - c_mfa(end))];
            end
            c_st_mfa = [c_st_mfa, c_mfa(end)];
        end

        % --- PA ---
        c_st_pa = [];
        for j = 1:length(f)
            c_pa = 0.0001;
            b_pa = 2*c_pa*(1 - c_pa);
            for t = 1:length(T)
                tu = b_pa(end)/(2*c_pa(end));
                td = b_pa(end)/(2*(1 - c_pa(end)));

                gp = (1 - c_pa(end))*((1 - p(pi))*td^q + p(pi)*f(j));
                gm = c_pa(end)*((1 - p(pi))*tu^q + p(pi)*(1 - f(j)));
                c_pa = [c_pa, c_pa(end) + dt*(gp - gm)];

                up = (2/k)*c_pa(end)*((1 - p(pi))*tu^q*(k - 2*q - 2*(k - q)*tu) + p(pi)*k*(1 - f(j))*(1 - 2*tu));
                down = (2/k)*(1 - c_pa(end))*((1 - p(pi))*td^q*(k - 2*q - 2*(k - q)*td) + p(pi)*k*f(j)*(1 - 2*td));
                b_pa = [b_pa, b_pa(end) + dt*(up + down)];
            end
            c_st_pa = [c_st_pa, c_pa(end)];
        end

        % --- SYMULACJE ---
        file_name = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_p_%d.txt', beta, k, q, round(p(pi)*100));
        if isfile(file_name)
            c_st_sim = csvread(file_name);
            c_st_sim_cut = c_st_sim(1:length(f));
        else
            warning('Brak pliku: %s', file_name);
            c_st_sim_cut = nan(1,length(f));
        end

        % dostosowanie próbkowania
        if k == 9999
            if p(pi) == 0.3
                c_sim_down = c_st_sim_cut(1:44);
                c_sim_middle = c_st_sim_cut(45:55);
                c_sim_up = c_st_sim_cut(56:end);
                c_st_sim_cut = [c_sim_down(1:2:end), c_sim_middle, c_sim_up(2:2:end)];
                f_sim = [f(1:2:44), f(45:55), f(56:2:end)];
            else
                c_st_sim_cut = c_st_sim_cut(1:2:end);
                f_sim = f(1:2:end);
            end
        elseif beta == 100
            c_st_sim_cut = c_st_sim_cut(1:2:end);
            f_sim = f(1:2:end);
        else
            if p(pi) == 0.1
                c_sim_down = c_st_sim_cut(1:44);
                c_sim_middle = c_st_sim_cut(45:55);
                c_sim_up = c_st_sim_cut(56:end);
                c_st_sim_cut = [c_sim_down(1:2:end), c_sim_middle, c_sim_up(2:2:end)];
                f_sim = [f(1:2:44), f(45:55), f(56:2:end)];
            else
                c_st_sim_cut = c_st_sim_cut(1:2:end);
                f_sim = f(1:2:end);
            end
        end

        % --- RYSOWANIE ---
        if points(pi) == 's'; marker_size = 6; else; marker_size = 5; end
        plot(f, c_st_pa, '-', 'Color', colors(pi,:), 'LineWidth', 1); hold on
        plot(f, c_st_mfa, '-', 'Color', colors(pi,:), 'LineWidth', 3); hold on
        LH(pi) = plot(f_sim, c_st_sim_cut, points(pi), 'MarkerFaceColor', colors(pi,:), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', marker_size); hold on

        leg{pi} = sprintf('$p^{ind} = %.1f$', p(pi));
    end

    % --- FORMATOWANIE ---
    ylim([0,1]); xlim([0,1]);
    axis square
    xlabel('$p_{eng}$', 'Interpreter','latex', 'FontSize', 20)
    ylabel('$c$', 'Interpreter','latex', 'FontSize', 20)
    h = legend(LH, leg, 'Location','northwest');
    set(h, 'Interpreter','latex', 'FontSize', 20)
end
