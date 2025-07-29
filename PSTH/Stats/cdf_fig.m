 function [ax, lgd] = cdf_fig(m1, m2, bno, l1, l2, xl)
%CDF_FIG Plot cumulative distribution functions (CDFs) for two datasets.
%
%   CDF_FIG(M1, M2, BNO, L1, L2, XL) plots the empirical cumulative 
%   distribution functions (CDFs) of two datasets M1 and M2 using histogram 
%   binning. The number of bins is specified by BNO. Legends are labeled as 
%   L1 and L2, and the x-axis label is set to XL.
%
%   Inputs:
%       m1 - First dataset (vector or array of numeric values)
%       m2 - Second dataset (vector or array of numeric values)
%       bno - Number of bins to use for histogram-based CDF estimation
%       l1 - Legend label for the first dataset (string)
%       l2 - Legend label for the second dataset (string)
%       xl - X-axis label and determines title; must match one of the following:
%            'Min FR', 'Max FR', 'Delay', 'Baseline', 'Cue', 'Reward',
%            'Punishment', 'Overall', 'First', 'Second'
%   Example:
%       cdf_fig(data1, data2, 50, 'Group A', 'Group B', 'Delay');
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Create shared bin edges based on min/max of both datasets
    edges = linspace(min(min(m1),min(m2)),max(max(m1),max(m2)),bno);
    
    figure;
    
    % Compute histogram counts for m1 and prepend zero to start from origin
    M1 = histcounts(m1, edges);
    M1= [0 M1(1:end)];  
    
    % Compute normalized cumulative sum for m1
    stairs(edges, cumsum(M1)./max(cumsum(M1)),'Color','r' , 'LineWidth', 1);
    hold on;
    
    % Same for M2
    M2 = histcounts(m2, edges);
    M2=[0 M2(1:end)];
    stairs(edges, cumsum(M2)./max(cumsum(M2)),'Color', 'b', 'LineWidth', 1);
    
    % Format figure
    axis tight
    lgd=legend(l1, l2);
    lgd.Location='best';
    lgd.FontSize=7;
    lgd.Box='off';
    xlabel(xl);
    
    if  strcmp('Min FR',xl)==1
        title({'CDF of'; 'Minimum FR (delay)'});
    elseif  strcmp('Max FR',xl)==1
        title({'CDF of'; 'Maximum FR (delay)'});
    elseif strcmp('Delay',xl)==1
        title({'CDF of'; 'Mean FR - Delay'});
    elseif  strcmp('Baseline',xl)==1
        title({'CDF of'; 'Mean FR - Baseline'});
    elseif strcmp('Cue',xl)==1
        title({'CDF of'; 'Mean FR - Cue'});
    elseif strcmp('Reward',xl)==1
        title({'CDF of'; 'Mean FR - Reward'});
    elseif strcmp('Punishment',xl)==1
        title({'CDF of'; 'Mean FR - Punishment'});
    elseif strcmp('Overall',xl)==1
        title({'CDF of'; 'Mean overall FR'});
    elseif strcmp('First',xl)==1
        title({'CDF of'; '1st Half - Delay'});
     elseif strcmp('Second',xl)==1
        title({'CDF of'; '2nd Half - Delay'});
    end
    
    ax=gca;
    ax.TickDir = 'out';
    ax.Box='off';
    hold off;

end