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
eco_stacj_cp_WS_analytical([0.25:0.25:0.75], 4, 9999, 100);

% --- SUBPLOT 2 ---
ax2 = axes('Position', [left_margins(2), 0.055, width, height]);
eco_stacj_cp_WS_analytical([0.25:0.25:0.75], 4, 8, 100);

% --- SUBPLOT 3 ---
ax3 = axes('Position', [left_margins(3), 0.055, width, height]);
eco_stacj_cp_WS_analytical([0.25:0.25:0.75], 4, 8, 5);

% --- ZAPIS DO PDF ---
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardcopy', 'off');
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches', ...
          'PaperSize', [pos(3) pos(4)], ...
          'PaperPosition', [0 0 pos(3) pos(4)]);
% print(fig, 'Cp_diff_p_eng_MFA_PA_MCS_arrows.pdf', '-dpdf', '-r300')


% =======================================================================
%  FUNKCJA: eco_stacj_cp_WS_analytical
% =======================================================================

function eco_stacj_cp_WS_analytical(p_eng, q, k, beta)
    % Funkcja rysuje stacjonarne wartoœci c(p)
    % dla podejœæ: MFA, PA oraz symulacji MC
    % na grafie Wattsa-Strogatza
    %
    % Argumenty:
    %   p_eng – wektor egnagement probability (np. [0.25 0.5 0.75])
    %   q     – rozmiar grupy wp³ywu (np. 4)
    %   k     – œrednia liczba s¹siadów (np. 8)
    %   beta  – rewiring probability w sieci Wattsa-Strogatza * 100
    %
    % Przyk³ad:
    %   eco_stacj_cp_WS_analytical([0.25:0.25:0.75], 4, 8, 100)

    points = ['^', 's', 'o', '*', 'd', 'h', 'v', 'x'];
    colors = [224, 159, 62; 158, 42, 43; 87, 117, 144]./255;

    LH = [];
    leg = [];
    for fi = 1:length(p_eng)
        if p_eng(fi) < 0.5
            % --- MFA ---
            c_mfa_up = 0.5:0.001:0.9999;
            p_st_mfa_up = (c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q)./ ...
                          (c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q + p_eng(fi) - c_mfa_up);

            c_mfa_down = 0.0001:0.001:0.4999;
            p_st_mfa_down = (c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q)./ ...
                            (c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q + p_eng(fi) - c_mfa_down);
            c_mfa_down = c_mfa_down(p_st_mfa_down >= 0 & p_st_mfa_down <= 1);
            p_st_mfa_down = p_st_mfa_down(p_st_mfa_down >= 0 & p_st_mfa_down <= 1);

            p_max_mfa = max(p_st_mfa_up);
            p_max_mfa_idx = find(p_st_mfa_up == p_max_mfa);

            c_mfa_unstable = c_mfa_up(1:p_max_mfa_idx);
            p_st_mfa_unstable = p_st_mfa_up(1:p_max_mfa_idx);
            c_mfa_stable = c_mfa_up(p_max_mfa_idx:end);
            p_st_mfa_stable = p_st_mfa_up(p_max_mfa_idx:end);

            % --- PA ---
            c_pa_up = 0.5:0.001:0.9999;
            b_st_pa_up = 2*(c_pa_up.*(1 - c_pa_up).*((1 - c_pa_up).^q * p_eng(fi) - c_pa_up.^q * (1 - p_eng(fi))) ...
                - (q/k)*(p_eng(fi) - c_pa_up).*(c_pa_up.*(1 - c_pa_up).^q + (1 - c_pa_up).*c_pa_up.^q)) ./ ...
                ((1 - c_pa_up).^q * p_eng(fi) - c_pa_up.^q * (1 - p_eng(fi)) - (q/k)*(p_eng(fi) - c_pa_up).*((1 - c_pa_up).^q + c_pa_up.^q));
            tu_up = b_st_pa_up./(2*c_pa_up);
            td_up = b_st_pa_up./(2*(1 - c_pa_up));
            p_st_pa_up = ((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q))./ ...
                         ((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q) - p_eng(fi) + c_pa_up);

            c_pa_down = 0.0001:0.001:0.4999;
            b_st_pa_down = 2*(c_pa_down.*(1 - c_pa_down).*((1 - c_pa_down).^q * p_eng(fi) - c_pa_down.^q * (1 - p_eng(fi))) ...
                - (q/k)*(p_eng(fi) - c_pa_down).*(c_pa_down.*(1 - c_pa_down).^q + (1 - c_pa_down).*c_pa_down.^q)) ./ ...
                ((1 - c_pa_down).^q * p_eng(fi) - c_pa_down.^q * (1 - p_eng(fi)) - (q/k)*(p_eng(fi) - c_pa_down).*((1 - c_pa_down).^q + c_pa_down.^q));
            tu_down = b_st_pa_down./(2*c_pa_down);
            td_down = b_st_pa_down./(2*(1 - c_pa_down));
            p_st_pa_down = ((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q))./ ...
                           ((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q) - p_eng(fi) + c_pa_down);
            c_pa_down = c_pa_down(p_st_pa_down >= 0 & p_st_pa_down <= 1);
            p_st_pa_down = p_st_pa_down(p_st_pa_down >= 0 & p_st_pa_down <= 1);

            p_max_pa = max(p_st_pa_up);
            p_max_pa_idx = find(p_st_pa_up == p_max_pa);
            c_pa_unstable = c_pa_up(1:p_max_pa_idx);
            p_st_pa_unstable = p_st_pa_up(1:p_max_pa_idx);
            c_pa_stable = c_pa_up(p_max_pa_idx:end);
            p_st_pa_stable = p_st_pa_up(p_max_pa_idx:end);

            % --- SYMULACJE ---
            p = 0:0.01:1;
            file_name = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_f_%d.txt', beta, k, q, p_eng(fi)*100);
            if isfile(file_name)
                c_st_sim = csvread(file_name);
                c_st_sim = c_st_sim(1:end-1);
                c_st_sim = c_st_sim(1:2:end);
                p_sim = p(1:2:end);
            else
                warning('Brak pliku: %s', file_name);
                c_st_sim = NaN(size(p(1:2:end)));
                p_sim = p(1:2:end);
            end

            % --- RYSOWANIE ---
            if points(fi) == 's'; marker_size = 6; else; marker_size = 5; end
            plot(p_st_pa_unstable, c_pa_unstable, '--', 'Color', colors(fi,:), 'LineWidth', 1); hold on
            plot(p_st_pa_stable, c_pa_stable, '-', 'Color', colors(fi,:), 'LineWidth', 1); hold on
            plot(p_st_pa_down, c_pa_down, '-', 'Color', colors(fi,:), 'LineWidth', 1); hold on
            plot(p_st_mfa_unstable, c_mfa_unstable, '--', 'Color', colors(fi,:), 'LineWidth', 3); hold on
            plot(p_st_mfa_stable, c_mfa_stable, '-', 'Color', colors(fi,:), 'LineWidth', 3); hold on
            plot(p_st_mfa_down, c_mfa_down, '-', 'Color', colors(fi,:), 'LineWidth', 3); hold on
            LH(fi) = plot(p_sim, c_st_sim, points(fi), 'MarkerFaceColor', colors(fi,:), 'MarkerEdgeColor', 'k', 'MarkerSize', marker_size); hold on
            leg{fi} = sprintf('$p_{eng} = %.2f$', p_eng(fi));
                elseif p_eng(fi) == 0.5
            % --- SYMULACJE ---
            p = 0:0.01:1;
            file_name = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_f_%d.txt', beta, k, q, p_eng(fi)*100);
            if isfile(file_name)
                c_st_sim = csvread(file_name);
                c_st_sim = c_st_sim(1:end-1);
                c_st_sim = c_st_sim(1:2:end);
                p_sim = p(1:2:end);
            else
                warning('Brak pliku: %s', file_name);
                c_st_sim = NaN(size(p(1:2:end)));
                p_sim = p(1:2:end);
            end

            % --- MFA i PA (symetryczne wzglêdem c = 0.5) ---
            c_mfa = 0:0.001:1;
            p_st_mfa = (c_mfa.*(1 - c_mfa).^q - (1 - c_mfa).*c_mfa.^q)./ ...
                        (c_mfa.*(1 - c_mfa).^q - (1 - c_mfa).*c_mfa.^q + 0.5 - c_mfa);

            c_pa = 0:0.001:1;
            b_st_pa = 2*(c_pa.*(1 - c_pa).*((1 - c_pa).^q * 0.5 - c_pa.^q * 0.5) ...
                - (q/k)*(0.5 - c_pa).*(c_pa.*(1 - c_pa).^q + (1 - c_pa).*c_pa.^q)) ./ ...
                ((1 - c_pa).^q * 0.5 - c_pa.^q * 0.5 - (q/k)*(0.5 - c_pa).*((1 - c_pa).^q + c_pa.^q));
            tu = b_st_pa./(2*c_pa);
            td = b_st_pa./(2*(1 - c_pa));
            p_st_pa = ((1 - c_pa).*(td.^q) - c_pa.*(tu.^q))./ ...
                      ((1 - c_pa).*(td.^q) - c_pa.*(tu.^q) - 0.5 + c_pa);

            % --- RYSOWANIE ---
            if points(fi) == 's'; marker_size = 6; else; marker_size = 5; end
            plot(p_st_pa, c_pa, '-', 'Color', colors(fi,:), 'LineWidth', 1); hold on
            plot(p_st_mfa, c_mfa, '-', 'Color', colors(fi,:), 'LineWidth', 3); hold on
            LH(fi) = plot(p_sim, c_st_sim, points(fi), 'MarkerFaceColor', colors(fi,:), ...
                          'MarkerEdgeColor', 'k', 'MarkerSize', marker_size); hold on
            leg{fi} = sprintf('$p_{eng} = %.2f$', p_eng(fi));

                elseif p_eng(fi) > 0.5
            % --- MFA ---
            c_mfa_up = 0.5001:0.001:0.9999;
            p_st_mfa_up = (c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q)./ ...
                          (c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q + p_eng(fi) - c_mfa_up);
            % filtrujemy dopuszczalne wartoœci p
            valid_up = (p_st_mfa_up >= 0 & p_st_mfa_up <= 1);
            c_mfa_up = c_mfa_up(valid_up);
            p_st_mfa_up = p_st_mfa_up(valid_up);

            c_mfa_down = 0.0001:0.001:0.5;
            p_st_mfa_down = (c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q)./ ...
                            (c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q + p_eng(fi) - c_mfa_down);

            % wyznaczamy "prze³¹cznik" (maksimum) na dolnej ga³êzi
            p_max_mfa = max(p_st_mfa_down);
            p_max_mfa_idx = find(p_st_mfa_down == p_max_mfa, 1, 'first');

            c_mfa_stable = c_mfa_down(1:p_max_mfa_idx);
            p_st_mfa_stable = p_st_mfa_down(1:p_max_mfa_idx);
            c_mfa_unstable = c_mfa_down(p_max_mfa_idx:end);
            p_st_mfa_unstable = p_st_mfa_down(p_max_mfa_idx:end);

            % --- PA ---
            c_pa_up = 0.5001:0.001:0.9999;
            b_st_pa_up = 2*(c_pa_up.*(1 - c_pa_up).*((1 - c_pa_up).^q * p_eng(fi) - c_pa_up.^q * (1 - p_eng(fi))) ...
                - (q/k)*(p_eng(fi) - c_pa_up).*(c_pa_up.*(1 - c_pa_up).^q + (1 - c_pa_up).*c_pa_up.^q)) ./ ...
                ((1 - c_pa_up).^q * p_eng(fi) - c_pa_up.^q * (1 - p_eng(fi)) - (q/k)*(p_eng(fi) - c_pa_up).*((1 - c_pa_up).^q + c_pa_up.^q));
            tu_up = b_st_pa_up./(2*c_pa_up);
            td_up = b_st_pa_up./(2*(1 - c_pa_up));
            p_st_pa_up = ((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q))./ ...
                         ((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q) - p_eng(fi) + c_pa_up);
            % filtrujemy dopuszczalne wartoœci p
            valid_pa_up = (p_st_pa_up >= 0 & p_st_pa_up <= 1);
            c_pa_up = c_pa_up(valid_pa_up);
            p_st_pa_up = p_st_pa_up(valid_pa_up);

            c_pa_down = 0.0001:0.001:0.5;
            b_st_pa_down = 2*(c_pa_down.*(1 - c_pa_down).*((1 - c_pa_down).^q * p_eng(fi) - c_pa_down.^q * (1 - p_eng(fi))) ...
                - (q/k)*(p_eng(fi) - c_pa_down).*(c_pa_down.*(1 - c_pa_down).^q + (1 - c_pa_down).*c_pa_down.^q)) ./ ...
                ((1 - c_pa_down).^q * p_eng(fi) - c_pa_down.^q * (1 - p_eng(fi)) - (q/k)*(p_eng(fi) - c_pa_down).*((1 - c_pa_down).^q + c_pa_down.^q));
            tu_down = b_st_pa_down./(2*c_pa_down);
            td_down = b_st_pa_down./(2*(1 - c_pa_down));
            p_st_pa_down = ((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q))./ ...
                           ((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q) - p_eng(fi) + c_pa_down);

            % wyznaczamy "prze³¹cznik" (maksimum) na dolnej ga³êzi PA
            p_max_pa = max(p_st_pa_down);
            p_max_pa_idx = find(p_st_pa_down == p_max_pa, 1, 'first');

            c_pa_stable = c_pa_down(1:p_max_pa_idx);
            p_st_pa_stable = p_st_pa_down(1:p_max_pa_idx);
            c_pa_unstable = c_pa_down(p_max_pa_idx:end);
            p_st_pa_unstable = p_st_pa_down(p_max_pa_idx:end);

            % --- SYMULACJE ---
            % zachowujemy specjalny przypadek dla p_eng==0.75 i beta==5 tak jak w oryginalnym kodzie
            if p_eng(fi) == 0.75 && beta == 5
                file_name_rzadko = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_f_%d.txt', beta, k, q, p_eng(fi)*100);
                if isfile(file_name_rzadko)
                    c_st_sim_rzadko = csvread(file_name_rzadko);
                    c_st_sim_rzadko = c_st_sim_rzadko(1:end-1);
                    c_st_sim_rzadko_up = c_st_sim_rzadko([0:0.01:1] > 0.05);
                else
                    c_st_sim_rzadko_up = NaN;
                end

                file_name_gesto = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_f_%d_gesto.txt', beta, k, q, p_eng(fi)*100);
                if isfile(file_name_gesto)
                    c_st_sim_gesto = csvread(file_name_gesto);
                    c_st_sim_gesto = c_st_sim_gesto(1:end-1);
                else
                    c_st_sim_gesto = NaN;
                end

                % ³¹czenie danych (tak jak w oryginale)
                if ~isnan(c_st_sim_gesto)
                    c_st_sim = [c_st_sim_gesto(1:14), c_st_sim_rzadko_up(1:2:end)];
                    p_rzadko = 0.06:0.01:1;
                    p_gesto = 0:0.002:0.05;
                    p = [p_gesto(1:14), p_rzadko(1:2:end)];
                else
                    p = 0:0.01:1;
                    file_name = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_f_%d.txt', beta, k, q, p_eng(fi)*100);
                    if isfile(file_name)
                        c_st_sim = csvread(file_name);
                        c_st_sim = c_st_sim(1:end-1);
                        c_st_sim = c_st_sim(1:2:end);
                    else
                        c_st_sim = NaN(size(0:0.01:1));
                    end
                end
            else
                p = 0:0.01:1;
                file_name = sprintf('stacjonarny_N_10000_beta_%d_k_%d_q_%d_f_%d.txt', beta, k, q, p_eng(fi)*100);
                if isfile(file_name)
                    c_st_sim = csvread(file_name);
                    c_st_sim = c_st_sim(1:end-1);
                    c_st_sim = c_st_sim(1:2:end);
                else
                    warning('Brak pliku: %s', file_name);
                    c_st_sim = NaN(size(p(1:2:end)));
                end
                p = p(1:2:end);
            end

            % --- RYSOWANIE ---
            if points(fi) == 's'; marker_size = 6; else; marker_size = 5; end
            % PA
            plot(p_st_pa_unstable, c_pa_unstable, '--', 'Color', colors(fi, :), 'LineWidth', 1); hold on
            plot(p_st_pa_stable, c_pa_stable, '-', 'Color', colors(fi, :), 'LineWidth', 1); hold on
            plot(p_st_pa_up, c_pa_up, '-', 'Color', colors(fi, :), 'LineWidth', 1); hold on
            % MFA
            plot(p_st_mfa_unstable, c_mfa_unstable, '--', 'Color', colors(fi, :), 'LineWidth', 3); hold on
            plot(p_st_mfa_stable, c_mfa_stable, '-', 'Color', colors(fi, :), 'LineWidth', 3); hold on
            plot(p_st_mfa_up, c_mfa_up, '-', 'Color', colors(fi, :), 'LineWidth', 3); hold on
            % SYMULACJE
            LH(fi) = plot(p, c_st_sim, points(fi), 'MarkerFaceColor', colors(fi, :), 'MarkerEdgeColor', 'k', 'MarkerSize', marker_size); hold on
            leg{fi} = sprintf('$p_{eng} = %0.2f$', p_eng(fi));
        end
    end

    ylim([0, 1]);
    xlim([0, 0.65]);
    xlabel('$$p^\mathrm{ind}$$', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('$$c$$', 'Interpreter', 'latex', 'FontSize', 20)
    axis square
    h = legend(LH, leg);
    set(h, 'Interpreter', 'latex', 'Location', 'east', 'FontSize', 20)
end
