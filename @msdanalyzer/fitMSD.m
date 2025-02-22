function obj = fitMSD(obj, clip_factor)
%%FITMSD Fit all MSD curves by a linear function.
%
% obj = obj.fitMSD fits all MSD curves by a straight line
%                      y = a * x + b.
% The fit is therefore rigorously valid only for purely
% diffusive behavior.
%
% Results are stored in the 'fit' field of the returned
% object. It is a structure with 2 fields:
% - a: all the values of the slope of the linear fit.
% - b: all the values for the intersect of the linear fit.
% - r2fit: the adjusted R2 value as a indicator of the goodness
% of the fit.
%
% obj = obj.fitMSD(clip_factor) does the fit, taking into
% account only the first potion of the average MSD curve
% specified by 'clip_factor' (a double between 0 and 1). If the
% value exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.


if nargin < 2
    clip_factor = 0.25;
end

if ~obj.msd_valid
    obj = obj.computeMSD;
end
n_spots = numel(obj.msd);    
    
msd_spot = obj.msd{1};

t = msd_spot(:,1);
y = msd_spot(:,2);
std = msd_spot(:,3);
w = 1./(std.^2);

% Clip data, never take the first one dt = 0
if clip_factor < 1
    t_limit = 2 : round(numel(t) * clip_factor);
else
    t_limit = 2 : min(1+round(clip_factor), numel(t));
end
t = t(t_limit);
y = y(t_limit);
w = w(t_limit);

% Thrash bad data
nonnan = ~isnan(y);
x = t(nonnan);
y = y(nonnan);
w = w(nonnan);

if numel(y) > 2
    % Compute sums needed for linear regression with errors
    Sxx = sum(w .* x.^2);
    Sxy = sum(w .* x .* y);
    
    % Compute best-fit parameters
    a = Sxy / Sxx; % Intercept
    
    % Compute uncertainties in parameters
    sigma_a = sqrt(1/Sxx);
        
    obj.lfit = struct(...
        'a', a, ...
        'err_a', sigma_a);
else 
    obj.lfit = struct(...
        'a', NaN, ...
        'err_a', NaN);
end
end