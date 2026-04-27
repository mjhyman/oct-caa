function mean_std_plot(x, exp_mu, exp_sd, ctl_mu, ctl_sd, ylims, tstr)
% mean_std_plot: Create a Mean ± Standard Deviation plot.
%   The STD will be a translucent ribbon surrounding the mean.
% 
% Inputs:
%   x (vector) - x-axis data points
%   exp_mu (vector) - vector of means for experimental
%   exp_sd (vector) - vector of standard deviations for experimental
%   ctl_mu (vector) - vector of means for control
%   ctl_sd (vector) - vector of standard deviations for control
%   ylims (string) - y-axis limits
%   tstr (string): title string

%% Esnure the x-axis and y-axis vectors are columns
% unwrap all into a column
x = x(:);
exp_mu = exp_mu(:);
exp_sd = exp_sd(:);
ctl_mu = ctl_mu(:);
ctl_sd = ctl_sd(:);

%% Set plotting parameters
% Set the experimental and control colors
exp_color = validatecolor('#DB5829');
ctl_color = validatecolor('#1964B0');
% Set the transparency level (0-1)
alpha = 0.3;

%% Create a figure and plot the mean
hold on;
% Scatter plots
scatter(x, exp_mu, 100, 'o', 'MarkerFaceColor',exp_color,...
        'DisplayName','Experimental Mean');
scatter(x, ctl_mu, 100, 'o', 'MarkerFaceColor',ctl_color,...
        'DisplayName', 'Control Mean');

%% Add ribbons
% Create ribbons for experimental and control means
ribbon(x, exp_mu, exp_sd, exp_color, alpha);
ribbon(x, ctl_mu, ctl_sd, ctl_color, alpha);

%% Add labels and legend
ylabel('Percentage Change');
ylim(ylims)
xlim([0,500]);
xlabel('Distance (\mum)')
% Add title and move higher
t = title(tstr);
pos = get(t,'Position');
set(t, 'Position', [pos(1), pos(2) + 1, pos(3)]);
hold off;
grid on;

%% Create Ribbon
    function ribbon(x, mu, sigma, color, alpha)      
        %%% Create x-axis points for ribbon
        x_lower = x;
        x_upper = flip(x);
        % Combine x_lower and x_upper then transpose
        x_fill = [x_lower; x_upper]';
        
        %%% Create y-axis points for ribbon
        y_lower = mu - sigma;
        y_upper = mu + sigma;        
        % Add extra points to y-xis to ensure vertical ends
        y_fill = [y_lower; flip(y_upper)]'; % y-coordinates for the fill
        
        %%% Fill +/-1 standard deviation b/w lower and upper
        fill(x_fill, y_fill, color,'FaceAlpha',alpha,'EdgeColor','none');        
    end
end
