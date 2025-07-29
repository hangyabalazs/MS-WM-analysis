function [mean_d1, se_d1, mean_d2, se_d2, Wp] = barmeanstat(d1, d2, l1, l2, alpha, str, parts)
% BARMEANSTAT   Bar plot with statistics.
%   BARMEANSTAT(D1, D2, L1, L2) plots two bar charts for the mean of
%   D1 and D2, with error bars indicating standard error. It performs
%   a Mann-Whitney U-test with significance level ALPHA (default is 0.05) 
%   and displays the significance level.
%
%   [MEAN_D1, SE_D1, MEAN_D2, SE_D2, WP] = BARMEANSTAT(D1, D2, L1, L2) returns
%   the mean and standard error for both datasets and the p-value from the
%   Mann-Whitney U-test.
%
%   BARMEANSTAT(___, ALPHA) specifies the significance level ALPHA for the
%   statistical test (default is 0.05).
%
%   BARMEANSTAT(___, STR) specifies the type of statistical test:
%       'nonpaired' - Mann-Whitney U-test (default)
%       'paired'    - Wilcoxon signed rank test
%
%   BARMEANSTAT(___, PARTS) controls the color scheme of the bar plot:
%       0 - Red (Exp) and Blue  (Ctrl)
%       1 - Green (Correct trials) and Red (Incorrect trials)
%
%   All input arguments can be used together in any combination that matches
%   the calling syntax above.
%
%   Example:
%       d1 = randn(100,1)+5;
%       d2 = randn(100,1)+4.8;
%       barmeanstat(d1, d2, 'Group A', 'Group B', 0.05, 'nonpaired', 1);
%
%   See also BOXPLOT.

    % Input argument check
    narginchk(4, 7);
    if nargin < 7
        str = 'nonpaired';
    end
    if nargin < 5 || isempty(alpha)
        alpha = 0.05; % default significance level
    end
    
    % Calculate means and standard errors
    mean_d1 = mean(d1);
    mean_d2 = mean(d2);
    se_d1 = std(d1) / sqrt(length(d1));
    se_d2 = std(d2) / sqrt(length(d2));
    
    % Bar plot
    figure;
    hold on;
    b = bar([1 2], [mean_d1, mean_d2], 'FaceColor', 'none', 'EdgeColor', 'flat', 'LineWidth', 1);
    if parts==1
      b.CData(1,:) = [0.3922    0.8314    0.0745]; % Green edge for d1
      b.CData(2,:) = [1 0 0]; % Red edge for d2  
    else
        b.CData(1,:) = [1 0 0]; % Red edge for d1
        b.CData(2,:) = [0 0 1]; % Blue edge for d2
    end
    
    % Add error bars
    errorbar([1 2], [mean_d1, mean_d2], [se_d1, se_d2], 'k', 'LineStyle', 'none', 'CapSize',0);
    
    % Perform statistical test
    switch str
        case 'nonpaired'
            [Wp, Wh] = ranksum(d1, d2, 'alpha', alpha); % Mann-Whitney U-test
        case 'paired'
            [Wp, Wh] = signrank(d1, d2, 'alpha', alpha); % Wilcoxon signed rank test
        otherwise
            error('boxstatmod:inputArg', 'Unsupported input argument.');
    end
    
    % Determine text color based on statistical significance
    clr = 'black';
    if Wh
        if Wp < 0.001
            star = '***';
        elseif Wp < 0.01
            star = '**';
        elseif Wp < alpha
            star = '*';
        end
    else
        star='ns';
    end
    
    % Add p-value text
    x_lim = xlim;
    tpos1 = (x_lim(1) + x_lim(2)) / 2;
    
    % Adjust y limits
    if mean_d1 >= 0 && mean_d2 >= 0
        m1=max(mean_d1 + se_d1, mean_d2 + se_d2);
        margin=m1*0.7;
        y_lim = [0, m1 + margin];
    elseif  mean_d1 <= 0 && mean_d2 <= 0
        m1=min(mean_d1 -se_d1, mean_d2 - se_d2);
        margin=m1*0.7;
        y_lim = [m1+margin, 0];
    elseif mean_d1 <= 0 && mean_d2 >= 0
        m1=min(mean_d1 -se_d1, mean_d2 - se_d2);
        m2=max(mean_d1 - se_d1, mean_d2 - se_d2);
        margin=m1*0.7;
        y_lim = [m1+margin, m2+margin];
    elseif mean_d1 >= 0 && mean_d2 <= 0
        m1=min(mean_d1 -se_d1, mean_d2 - se_d2);
        m2=max(mean_d1 - se_d1, mean_d2 - se_d2);
        margin=m1*0.7;
        y_lim = [m1+margin, m2+margin];
    end
    ylim(y_lim);
    
    % Plot significance stars
    if mean_d1 >= 0 && mean_d2 >= 0
        tpos2 = m1 + (m1*0.2);
        text(tpos1, tpos2, star, 'Color', clr, 'HorizontalAlignment', 'center', 'FontSize', 12);
        v1=m1+(abs(m1)*0.016);
        line([x_lim(1) + 1.2, x_lim(1) + 1.2], [tpos2 - (abs(m1)*0.0083), v1], 'Color', 'black');
        line([x_lim(1)+2.2,x_lim(1)+2.2], [tpos2 - (abs(m1)*0.0083), v1], 'Color', 'black');
        v=tpos2 - (abs(m1)*0.0083);
        line([x_lim(1) + 1.2, x_lim(1) + 1.2, x_lim(1)+2.2, x_lim(1)+2.2], [v,v,v,v], 'Color', 'black');
    else
        tpos2 = m1 + (m1*0.2) ;
        text(tpos1, tpos2 - (abs(m1)*0.0833), star, 'Color', clr, 'HorizontalAlignment', 'center', 'FontSize', 12);
        v1=m1-(abs(m1)*0.0775);    
        line([x_lim(1) + 2.2, x_lim(1) + 2.2], [tpos2 - (abs(m1)*0.0083), v1], 'Color', 'black');
        line([x_lim(1)+1.2,x_lim(1)+1.2], [tpos2 - (abs(m1)*0.0083), v1], 'Color', 'black');
        v=tpos2 - (abs(m1)*0.0083);
        line([x_lim(1) + 1.2, x_lim(1) + 1.2, x_lim(1)+2.2, x_lim(1)+2.2], [v,v,v,v], 'Color', 'black');
    end
    
    
    ax = gca;
    
    % Add custom legend using patches
    if parts == 1
        p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'green', 'LineWidth', 1);
        p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'red', 'LineWidth', 1);
    else 
        p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'red', 'LineWidth', 1);
        p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'blue', 'LineWidth', 1);
    end
    legend([p1, p2], {l1, l2}, 'Location', 'best','Box','off');
    
    % Format plot
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.XTick = [1 2];
    ax.XTickLabel = {l1, l2};
    hold off;

end
