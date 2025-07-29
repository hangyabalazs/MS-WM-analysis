function p_str = format_p_value(p)
%FORMAT_P_VALUE Convert p-value to significance stars
%
% This function converts a p-value to asterisk notation for statistical
% significance levels. The conversion follows standard conventions:
%   p < 0.001  ->  '***'
%   p < 0.01   ->  '**'
%   p < 0.05   ->  '*'
%   p >= 0.05  ->  'ns'
%
% Inputs:
%   p - scalar numeric p-value (should be between 0 and 1)
%
% Outputs:
%   p_str - character vector containing significance stars or empty string
%
% Examples:
%   format_p_value(0.0005)  % Returns '***'
%   format_p_value(0.005)   % Returns '**'
%   format_p_value(0.03)    % Returns '*'
%   format_p_value(0.1)     % Returns 'ns'

threestars=0.001;
twostars=0.01;

    if p < threestars
         p_str= '***';
    elseif p < twostars
            p_str= '**';
    elseif  p < 0.05
            p_str = '*';
    else 
        p_str ='ns';
    end
end
