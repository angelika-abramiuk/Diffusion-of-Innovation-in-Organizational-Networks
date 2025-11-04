% PARAMETRY OGÓLNE
colors = [224, 159, 62;
          158, 42, 43;
          21, 109, 160;
          34, 140, 34;
          128, 93, 179;
          13, 170, 187] ./ 255;

q_list = [2, 4, 6, 8];
p_eng_values = [0.55, 0.75];
legend_pos = {'southeast', 'southeast', 'southeast'};

% TWORZENIE FIGURY
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

% PANEL 1 i 2: MFA dla ró¿nych q i p_eng
for j = 1:2
    ax = axes('Position', [left_margins(j), 0.055, width, height]);
    hold(ax, 'on');
    box(ax, 'on');

    p_eng = p_eng_values(j);
    LH = gobjects(1, numel(q_list));
    leg = cell(1, numel(q_list));

    for qi = length(q_list):-1:1
        q = q_list(qi);

        if p_eng < 0.5
            % MFA (p_eng < 0.5) - zgodnie z Twoim orygina³em
            c_mfa_down = 0.0001:0.0001:0.4999;
            num = (c_mfa_down .* (1 - c_mfa_down).^q - (1 - c_mfa_down) .* c_mfa_down.^q);
            den = (c_mfa_down .* (1 - c_mfa_down).^q - (1 - c_mfa_down) .* c_mfa_down.^q + p_eng - c_mfa_down);
            p_st_mfa_down = num ./ den;

            % usuñ NaN i wartoœci spoza [0, 1]
            valid = isfinite(p_st_mfa_down) & (p_st_mfa_down >= 0) & (p_st_mfa_down <= 1);
            p_st_mfa_down = p_st_mfa_down(valid);
            c_mfa_down = c_mfa_down(valid);

            % znajdŸ lokalne maksima i minima
            [pks_max, locs_max] = findpeaks(p_st_mfa_down);
            [pks_min, locs_min] = findpeaks(-p_st_mfa_down);
            pks_min = -pks_min;

            % górna ga³¹Ÿ (rysowana zawsze)
            c_mfa_up = 0.5:0.0001:0.9999;
            p_st_mfa_up = (c_mfa_up .* (1 - c_mfa_up).^q - (1 - c_mfa_up) .* c_mfa_up.^q) ...
                           ./ (c_mfa_up .* (1 - c_mfa_up).^q - (1 - c_mfa_up) .* c_mfa_up.^q + p_eng - c_mfa_up);

            % stabilna / niestabilna czêœæ górna
            p_max_mfa_up = max(p_st_mfa_up);
            p_max_mfa_up_idx = find(p_st_mfa_up == p_max_mfa_up, 1);
            c_mfa_stable_up = c_mfa_up(p_max_mfa_up_idx:end);
            p_st_mfa_stable_up = p_st_mfa_up(p_max_mfa_up_idx:end);
            c_mfa_unstable_up = c_mfa_up(1:p_max_mfa_up_idx);
            p_st_mfa_unstable_up = p_st_mfa_up(1:p_max_mfa_up_idx);

            % rysowanie górnej ga³êzi
            plot(ax, p_st_mfa_unstable_up, c_mfa_unstable_up, '--', 'Color', colors(qi,:), 'LineWidth', 2);
            LH(qi) = plot(ax, p_st_mfa_stable_up, c_mfa_stable_up, '-', 'Color', colors(qi,:), 'LineWidth', 2);

            % dolna ga³¹Ÿ: warunek histerezy
            if isempty(pks_max) || isempty(pks_min)
                % brak histerezy - rysujemy ca³oœæ ci¹g³¹ lini¹
                plot(ax, p_st_mfa_down, c_mfa_down, '-', 'Color', colors(qi,:), 'LineWidth', 2);
            else
                % histereza - wydzielamy stabilne i niestabilne fragmenty
                [~, imax] = max(pks_max);
                [~, imin] = min(pks_min);
                c_max = c_mfa_down(locs_max(imax));
                c_min = c_mfa_down(locs_min(imin));

                i1 = find(c_mfa_down >= min(c_min, c_max), 1, 'first');
                i2 = find(c_mfa_down <= max(c_min, c_max), 1, 'last');

                % podzia³ na fragmenty
                c_mfa_stable_down_up = c_mfa_down(i2:end);
                p_st_mfa_stable_down_up = p_st_mfa_down(i2:end);

                c_mfa_stable_down_down = c_mfa_down(1:i1);
                p_st_mfa_stable_down_down = p_st_mfa_down(1:i1);

                c_mfa_unstable_down = c_mfa_down(i1:i2);
                p_st_mfa_unstable_down = p_st_mfa_down(i1:i2);

                % rysowanie dolnej czêœci
                plot(ax, p_st_mfa_unstable_down, c_mfa_unstable_down, '--', 'Color', colors(qi,:), 'LineWidth', 2);
                plot(ax, p_st_mfa_stable_down_up, c_mfa_stable_down_up, '-', 'Color', colors(qi,:), 'LineWidth', 2);
                plot(ax, p_st_mfa_stable_down_down, c_mfa_stable_down_down, '-', 'Color', colors(qi,:), 'LineWidth', 2);
            end

        elseif p_eng > 0.5
            % MFA (p_eng > 0.5)
            c_mfa_up = 0.5001:0.0001:0.9999;
            num = (c_mfa_up .* (1 - c_mfa_up).^q - (1 - c_mfa_up) .* c_mfa_up.^q);
            den = (c_mfa_up .* (1 - c_mfa_up).^q - (1 - c_mfa_up) .* c_mfa_up.^q + p_eng - c_mfa_up);
            p_st_mfa_up = num ./ den;

            valid = isfinite(p_st_mfa_up) & (p_st_mfa_up >= 0) & (p_st_mfa_up <= 1);
            p_st_mfa_up = p_st_mfa_up(valid);
            c_mfa_up = c_mfa_up(valid);

            [pks_max, locs_max] = findpeaks(p_st_mfa_up);
            [pks_min, locs_min] = findpeaks(-p_st_mfa_up);
            pks_min = -pks_min;

            % dolna ga³¹Ÿ
            c_mfa_down = 0.0001:0.0001:0.5;
            p_st_mfa_down = (c_mfa_down .* (1 - c_mfa_down).^q - (1 - c_mfa_down) .* c_mfa_down.^q) ...
                ./ (c_mfa_down .* (1 - c_mfa_down).^q - (1 - c_mfa_down) .* c_mfa_down.^q + p_eng - c_mfa_down);

            p_max_mfa_down = max(p_st_mfa_down);
            p_max_mfa_down_idx = find(p_st_mfa_down == p_max_mfa_down, 1);
            c_mfa_stable_down = c_mfa_down(1:p_max_mfa_down_idx);
            p_st_mfa_stable_down = p_st_mfa_down(1:p_max_mfa_down_idx);
            c_mfa_unstable_down = c_mfa_down(p_max_mfa_down_idx:end);
            p_st_mfa_unstable_down = p_st_mfa_down(p_max_mfa_down_idx:end);

            plot(ax, p_st_mfa_unstable_down, c_mfa_unstable_down, '--', 'Color', colors(qi,:), 'LineWidth', 2);
            LH(qi) = plot(ax, p_st_mfa_stable_down, c_mfa_stable_down, '-', 'Color', colors(qi,:), 'LineWidth', 2);

            % górna ga³¹Ÿ histerezy
            if isempty(pks_max) || isempty(pks_min)
                plot(ax, p_st_mfa_up, c_mfa_up, '-', 'Color', colors(qi,:), 'LineWidth', 2);
            else
                [~, imax] = max(pks_max);
                [~, imin] = min(pks_min);
                c_max = c_mfa_up(locs_max(imax));
                c_min = c_mfa_up(locs_min(imin));
                i1 = find(c_mfa_up >= min(c_min, c_max), 1, 'first');
                i2 = find(c_mfa_up <= max(c_min, c_max), 1, 'last');
                c_mfa_unstable_up = c_mfa_up(i1:i2);
                p_st_mfa_unstable_up = p_st_mfa_up(i1:i2);
                c_mfa_stable_up_up = c_mfa_up(i2:end);
                p_st_mfa_stable_up_up = p_st_mfa_up(i2:end);
                c_mfa_stable_up_down = c_mfa_up(1:i1);
                p_st_mfa_stable_up_down = p_st_mfa_up(1:i1);

                plot(ax, p_st_mfa_unstable_up, c_mfa_unstable_up, '--', 'Color', colors(qi,:), 'LineWidth', 2);
                plot(ax, p_st_mfa_stable_up_up, c_mfa_stable_up_up, '-', 'Color', colors(qi,:), 'LineWidth', 2);
                plot(ax, p_st_mfa_stable_up_down, c_mfa_stable_up_down, '-', 'Color', colors(qi,:), 'LineWidth', 2);
            end

        else
            % MFA dla p_eng = 0.5
            c_mfa = 0.0001:0.0001:0.9999;
            p_st_mfa = (c_mfa.*(1 - c_mfa).^q - (1 - c_mfa).*c_mfa.^q) ...
                      ./ (c_mfa.*(1 - c_mfa).^q - (1 - c_mfa).*c_mfa.^q + p_eng - c_mfa);

            % filtrujemy punkty w zakresie [0, 1]
            valid = p_st_mfa >= 0 & p_st_mfa <= 1;
            c_mfa = c_mfa(valid);
            p_st_mfa = p_st_mfa(valid);

            % szukamy ekstremów
            dpdc = gradient(p_st_mfa, c_mfa);
            extrema_idx = find(diff(sign(dpdc)));

            if numel(extrema_idx) >= 2
                % przypadek z histerez¹
                i1 = extrema_idx(1);
                i2 = extrema_idx(end);

                % ga³¹Ÿ dolna stabilna
                LH(qi) = plot(ax, p_st_mfa(1:i1), c_mfa(1:i1), '-', ...
                              'Color', colors(qi, :), 'LineWidth', 2);
                % ga³¹Ÿ niestabilna
                plot(ax, p_st_mfa(i1:i2), c_mfa(i1:i2), '--', ...
                     'Color', colors(qi, :), 'LineWidth', 2);
                % ga³¹Ÿ górna stabilna
                plot(ax, p_st_mfa(i2:end), c_mfa(i2:end), '-', ...
                     'Color', colors(qi, :), 'LineWidth', 2);

                % linia w c = 0.5
                p_min = min(p_st_mfa(i1:i2));
                plot(ax, [0 p_min], [0.5 0.5], '--', ...
                     'Color', colors(qi, :), 'LineWidth', 2);
                plot(ax, [p_min 1], [0.5 0.5], '-', ...
                     'Color', colors(qi, :), 'LineWidth', 2);

            else
                % brak histerezy
                LH(qi) = plot(ax, p_st_mfa, c_mfa, '-', ...
                              'Color', colors(qi, :), 'LineWidth', 2);

                % linia w c = 0.5
                p_max = max(p_st_mfa);
                plot(ax, [0 p_max], [0.5 0.5], '--', ...
                     'Color', colors(qi, :), 'LineWidth', 2);
                plot(ax, [p_max 1], [0.5 0.5], '-', ...
                     'Color', colors(qi, :), 'LineWidth', 2);
            end
        end

        leg{qi} = sprintf('$q = %0.1f$', q);
    end

    xlabel(ax, '$p^{\mathrm{ind}}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(ax, '$c$', 'Interpreter', 'latex', 'FontSize', 20);
    axis(ax, 'square');
    ylim(ax, [0, 1]);
    xlim(ax, [0, 0.45]);
    legend(ax, LH, leg, 'Interpreter', 'latex', 'FontSize', 18, 'Location', legend_pos{j});
end

% PANEL 3: PA i MFA dla ró¿nych k
ax3 = axes('Position', [left_margins(3), 0.055, width, height]);
hold(ax3, 'on');
box(ax3, 'on');

k = [8, 10, 20, 60, 100, 150];
q = 4;
p_eng = 0.75;

LH = [];
leg = [];
for ki = 1:length(k)
    if p_eng < 0.5
        % MFA
        c_mfa_up = [0.5:0.001:0.9999];
        p_st_mfa_up = (c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q)./(c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q + p_eng - c_mfa_up);

        c_mfa_down = [0.0001:0.001:0.4999];
        p_st_mfa_down = (c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q)./(c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q + p_eng - c_mfa_down);
        c_mfa_down = c_mfa_down(p_st_mfa_down >= 0 & p_st_mfa_down <= 1);
        p_st_mfa_down = p_st_mfa_down(p_st_mfa_down >= 0 & p_st_mfa_down <= 1);

        p_max_mfa = max(p_st_mfa_up);
        p_max_mfa_idx = find(p_st_mfa_up == p_max_mfa);

        c_mfa_unstable = c_mfa_up(1:p_max_mfa_idx);
        p_st_mfa_unstable = p_st_mfa_up(1:p_max_mfa_idx);

        c_mfa_stable = c_mfa_up(p_max_mfa_idx:end);
        p_st_mfa_stable = p_st_mfa_up(p_max_mfa_idx:end);

        % PA
        c_pa_up = [0.5:0.001:0.9999];
        b_st_pa_up = 2*(c_pa_up.*(1 - c_pa_up).*((1 - c_pa_up).^q * p_eng - c_pa_up.^q * (1 - p_eng)) - (q/k(ki))*(p_eng - c_pa_up).*(c_pa_up.*(1 - c_pa_up).^q + (1 - c_pa_up).*c_pa_up.^q))./((1 - c_pa_up).^q * p_eng - c_pa_up.^q * (1 - p_eng) - (q/k(ki))*(p_eng - c_pa_up).*((1 - c_pa_up).^q + c_pa_up.^q));
        tu_up = b_st_pa_up./(2*c_pa_up);
        td_up = b_st_pa_up./(2*(1 - c_pa_up));
        p_st_pa_up = ((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q))./((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q) - p_eng + c_pa_up);

        c_pa_down = [0.0001:0.001:0.4999];
        b_st_pa_down = 2*(c_pa_down.*(1 - c_pa_down).*((1 - c_pa_down).^q * p_eng - c_pa_down.^q * (1 - p_eng)) - (q/k(ki))*(p_eng - c_pa_down).*(c_pa_down.*(1 - c_pa_down).^q + (1 - c_pa_down).*c_pa_down.^q))./((1 - c_pa_down).^q * p_eng - c_pa_down.^q * (1 - p_eng) - (q/k(ki))*(p_eng - c_pa_down).*((1 - c_pa_down).^q + c_pa_down.^q));
        tu_down = b_st_pa_down./(2*c_pa_down);
        td_down = b_st_pa_down./(2*(1 - c_pa_down));
        p_st_pa_down = ((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q))./((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q) - p_eng + c_pa_down);
        c_pa_down = c_pa_down(p_st_pa_down >= 0 & p_st_pa_down <= 1);
        p_st_pa_down = p_st_pa_down(p_st_pa_down >= 0 & p_st_pa_down <= 1);

        p_max_pa = max(p_st_pa_up);
        p_max_pa_idx = find(p_st_pa_up == p_max_pa);

        c_pa_unstable = c_pa_up(1:p_max_pa_idx);
        p_st_pa_unstable = p_st_pa_up(1:p_max_pa_idx);

        c_pa_stable = c_pa_up(p_max_pa_idx:end);
        p_st_pa_stable = p_st_pa_up(p_max_pa_idx:end);

        % RYSOWANIE
        % PA
        plot(p_st_pa_unstable, c_pa_unstable, '--', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        plot(p_st_pa_stable, c_pa_stable, '-', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        LH(ki) = plot(p_st_pa_down, c_pa_down, '-', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        leg{ki} = sprintf('$$\\langle k \\rangle = %.0f$$', k(ki));
        % MFA
        plot(p_st_mfa_unstable, c_mfa_unstable, '--', 'Color', 'k', 'LineWidth', 2);
        hold on
        plot(p_st_mfa_stable, c_mfa_stable, '-', 'Color', 'k', 'LineWidth', 2);
        hold on
        LH(length(k) + 1) = plot(p_st_mfa_down, c_mfa_down, '-', 'Color', 'k', 'LineWidth', 2);
        hold on

    elseif f > 0.5
        % MFA
        c_mfa_up = [0.5001:0.001:0.9999];
        p_st_mfa_up = (c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q)./(c_mfa_up.*(1 - c_mfa_up).^q - (1 - c_mfa_up).*c_mfa_up.^q + p_eng - c_mfa_up);
        c_mfa_up = c_mfa_up(p_st_mfa_up >= 0 & p_st_mfa_up <= 1);
        p_st_mfa_up = p_st_mfa_up(p_st_mfa_up >= 0 & p_st_mfa_up <= 1);

        c_mfa_down = [0.0001:0.001:0.5];
        p_st_mfa_down = (c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q)./(c_mfa_down.*(1 - c_mfa_down).^q - (1 - c_mfa_down).*c_mfa_down.^q + p_eng - c_mfa_down);

        p_max_mfa = max(p_st_mfa_down);
        p_max_mfa_idx = find(p_st_mfa_down == p_max_mfa);

        c_mfa_stable = c_mfa_down(1:p_max_mfa_idx);
        p_st_mfa_stable = p_st_mfa_down(1:p_max_mfa_idx);

        c_mfa_unstable = c_mfa_down(p_max_mfa_idx:end);
        p_st_mfa_unstable = p_st_mfa_down(p_max_mfa_idx:end);

        % PA
        c_pa_up = [0.5001:0.001:0.9999];
        b_st_pa_up = 2*(c_pa_up.*(1 - c_pa_up).*((1 - c_pa_up).^q * p_eng - c_pa_up.^q * (1 - p_eng)) - (q/k(ki))*(p_eng - c_pa_up).*(c_pa_up.*(1 - c_pa_up).^q + (1 - c_pa_up).*c_pa_up.^q))./((1 - c_pa_up).^q * p_eng - c_pa_up.^q * (1 - p_eng) - (q/k(ki))*(p_eng - c_pa_up).*((1 - c_pa_up).^q + c_pa_up.^q));
        tu_up = b_st_pa_up./(2*c_pa_up);
        td_up = b_st_pa_up./(2*(1 - c_pa_up));
        p_st_pa_up = ((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q))./((1 - c_pa_up).*(td_up.^q) - c_pa_up.*(tu_up.^q) - p_eng + c_pa_up);
        c_pa_up = c_pa_up(p_st_pa_up >= 0 & p_st_pa_up <= 1);
        p_st_pa_up = p_st_pa_up(p_st_pa_up >= 0 & p_st_pa_up <= 1);

        c_pa_down = [0.0001:0.001:0.5];
        b_st_pa_down = 2*(c_pa_down.*(1 - c_pa_down).*((1 - c_pa_down).^q * p_eng - c_pa_down.^q * (1 - p_eng)) - (q/k(ki))*(p_eng - c_pa_down).*(c_pa_down.*(1 - c_pa_down).^q + (1 - c_pa_down).*c_pa_down.^q))./((1 - c_pa_down).^q * p_eng - c_pa_down.^q * (1 - p_eng) - (q/k(ki))*(p_eng - c_pa_down).*((1 - c_pa_down).^q + c_pa_down.^q));
        tu_down = b_st_pa_down./(2*c_pa_down);
        td_down = b_st_pa_down./(2*(1 - c_pa_down));
        p_st_pa_down = ((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q))./((1 - c_pa_down).*(td_down.^q) - c_pa_down.*(tu_down.^q) - p_eng + c_pa_down);

        p_max_pa = max(p_st_pa_down);
        p_max_pa_idx = find(p_st_pa_down == p_max_pa);

        c_pa_stable = c_pa_down(1:p_max_pa_idx);
        p_st_pa_stable = p_st_pa_down(1:p_max_pa_idx);

        c_pa_unstable = c_pa_down(p_max_pa_idx:end);
        p_st_pa_unstable = p_st_pa_down(p_max_pa_idx:end);

        % RYSOWANIE
        % PA
        plot(p_st_pa_unstable, c_pa_unstable, '--', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        plot(p_st_pa_stable, c_pa_stable, '-', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        LH(ki) = plot(p_st_pa_up, c_pa_up, '-', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        leg{ki} = sprintf('$$\\langle k \\rangle = %.0f$$', k(ki));
        % MFA
        plot(p_st_mfa_unstable, c_mfa_unstable, '--', 'Color', 'k', 'LineWidth', 2);
        hold on
        plot(p_st_mfa_stable, c_mfa_stable, '-', 'Color', 'k', 'LineWidth', 2);
        hold on
        LH(length(k) + 1) = plot(p_st_mfa_up, c_mfa_up, '-', 'Color', 'k', 'LineWidth', 2);
        hold on

    else
        % MFA
        c_mfa = [0.0001:0.001:0.9999];
        p_st_mfa = (c_mfa.*(1 - c_mfa).^q - (1 - c_mfa).*c_mfa.^q)./(c_mfa.*(1 - c_mfa).^q - (1 - c_mfa).*c_mfa.^q + p_eng - c_mfa);
        p_max_mfa = max(p_st_mfa);
        p_max_mfa_idx = find(p_st_mfa == p_max_mfa);

        % PA
        c_pa = [0.0001:0.001:0.9999];
        b_st_pa = 2*(c_pa.*(1 - c_pa).*((1 - c_pa).^q * f - c_pa.^q * (1 - p_eng)) - (q/k(ki))*(p_eng - c_pa).*(c_pa.*(1 - c_pa).^q + (1 - c_pa).*c_pa.^q))./((1 - c_pa).^q * p_eng - c_pa.^q * (1 - p_eng) - (q/k(ki))*(p_eng - c_pa).*((1 - c_pa).^q + c_pa.^q));
        tu = b_st_pa./(2*c_pa);
        td = b_st_pa./(2*(1 - c_pa));
        p_st_pa = ((1 - c_pa).*(td.^q) - c_pa.*(tu.^q))./((1 - c_pa).*(td.^q) - c_pa.*(tu.^q) - p_eng + c_pa);
        p_max_pa = max(p_st_pa);
        p_max_pa_idx = find(p_st_pa == p_max_pa);

        % RYSOWANIE
        % PA
        LH(ki) = plot(p_st_pa, c_pa, '-', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        plot([0, p_max_pa], [0.5, 0.5], '--', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        plot([p_max_pa, 1], [0.5, 0.5], '-', 'Color', colors(ki, :), 'LineWidth', 2);
        hold on
        leg{ki} = sprintf('$$\\langle k \\rangle = %.0f$$', k(ki));
        % MFA
        LH(length(k) + 1) = plot(p_st_mfa, c_mfa, '-', 'Color', 'k', 'LineWidth', 2);
        hold on
        plot([0, p_max_mfa], [0.5, 0.5], '--', 'Color', 'k', 'LineWidth', 2);
        hold on
        plot([p_max_mfa, 1], [0.5, 0.5], '-', 'Color', 'k', 'LineWidth', 2);
        hold on
    end
end

leg{ki + 1} = 'MFA';

ylim([0, 1]);
xlim([0, 0.4]);
xlabel('$$p^\mathrm{ind}$$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$$c$$', 'Interpreter', 'latex', 'FontSize', 20)

axis square

h = legend(LH, leg);
set(h, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 20)

% ZAPIS DO PDF (odkomentuj jeœli chcesz zapisaæ)
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardcopy', 'off');
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
set(fig, 'PaperUnits', 'inches', ...
          'PaperSize', [pos(3) pos(4)], ...
          'PaperPosition', [0 0 pos(3) pos(4)]);
% print(fig, 'MFA_PA_three_panel_fixed.pdf', '-dpdf', '-r300');
