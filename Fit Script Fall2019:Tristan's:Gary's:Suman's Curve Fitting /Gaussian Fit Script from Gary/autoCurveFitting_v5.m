% autoCurveFitting_v5.m

% This program will automatically fit both the fall and rise sections of
% the data after cutting is done.  The amplitudes and their uncertainties
% are saved as well as the exponents.

[filename, pathname, filterindex] = uigetfile({'*.trc',  'Matlab data files (*.mat)'; '*.*', 'All Files (*.*)'}, 'Select a data set', 'MultiSelect', 'on'); % load files
N = length(filename); % Number of files loaded

%Initialize Arrays for Extracted Exponentials
fall_smallexp = zeros(N,1); % Is the slow fall rate coefficients
fall_mediumexp = zeros(N,1); % Is the medium fall rate coefficients
fall_largeexp = zeros(N,1); % Is the fast fall rate coefficients
rise_smallexp = zeros(N,1); % Is the slow rise rate coefficients
rise_largeexp = zeros(N,1); % Is the fast rise rate coefficients

fall_smallexp_lowbnd = zeros(N,1); % Is the lower bound for slow fall rate coefficients
fall_mediumexp_lowbnd = zeros(N,1); % Is the lower bound for medium fall rate coefficients
fall_largeexp_lowbnd = zeros(N,1); % Is the lower bound for large fall rate coefficients
rise_smallexp_lowbnd = zeros(N,1); % Is the lowerbound for the small rise rate coefficients
rise_largeexp_lowbnd = zeros(N,1); % Is the lowerbound for the large rise rate coefficients

fall_smallexp_error = zeros(N,1); % Is the error for the slow fall rate coefficients
fall_mediumexp_error = zeros(N,1); % Is the error for the medium fall rate coefficients
fall_largeexp_error = zeros(N,1); % Is the error for the fast fall rate coefficients
rise_smallexp_error = zeros(N,1); % Is the error for the slow rise rate coefficients
rise_largeexp_error = zeros(N,1); % Is the error for the fast rise rate coefficients

%Initialize Arrays for Extracted Amplitudes
fall_smallamp = zeros(N,1); % Is the slow fall rate amplitude
fall_mediumamp = zeros(N,1); % Is the medium fall rate amplitude
fall_largeamp = zeros(N,1); % Is the fast fall rate amplitude
rise_smallamp = zeros(N,1); % Is the slow rise rate amplitude
rise_largeamp = zeros(N,1); % Is the fast rise rate amplitude

fall_smallamp_lowbnd = zeros(N,1); % Is the lower bound for slow fall rate amplitude
fall_mediumamp_lowbnd = zeros(N,1); % Is the lower bound for medium fall rate amplitude
fall_largeamp_lowbnd = zeros(N,1); % Is the lower bound for large fall rate amplitude
rise_smallamp_lowbnd = zeros(N,1); % Is the lowerbound for the small rise rate amplitude
rise_largeamp_lowbnd = zeros(N,1); % Is the lowerbound for the large rise rate amplitude

fall_smallamp_error = zeros(N,1); % Is the error for the slow fall rate amplitude
fall_mediumamp_error = zeros(N,1); % Is the error for the medium fall rate amplitude
fall_largeamp_error = zeros(N,1); % Is the error for the fast fall rate amplitude
rise_smallamp_error = zeros(N,1); % Is the error for the slow rise rate amplitude
rise_largeamp_error = zeros(N,1); % Is the error for the fast rise rate amplitude

fitresultarray = cell(N,3); % Contains all information for the fits

% Establish save path name
savepath = '\\ecas.wvu.edu\squol\Data\Charging and Relaxation Dynamics (Gary)\Trion\Data and Analysis\HeNe Modulation\Third Analysis';

for i=1:N
    % Load data
    load([pathname,filename{i}]);
    
    % Fit current data
    [fitresult, gof, hf1, hf2] = createFits_v5(tfall, yfall, trise, yrise);   % The main function to call to fit the data
    
    % Extract title for plot and save figures
    temp_HeNeOD = filename{i}(cumsum(filename{i}(:) == '_') == 1);
    HeNeOD = temp_HeNeOD(2:11);
    set(hf1,'Name',HeNeOD)
    set(hf2,'Name',HeNeOD)
    saveas(hf1,[savepath,HeNeOD,'_Fall'],'fig');
    saveas(hf2,[savepath,HeNeOD,'_Rise'],'fig');
    
    % Extract CTL Power to include in the name of the final file to save
    filename_string = filename{i}(cumsum(filename{i}(:) == '_') == 0);
    filename_prefix = filename_string(16:22);
    
    %Establish arrays to contain all fit parameters
    fitresultarray{i,1} = fitresult{1};
    fitresultarray{i,2} = fitresult{2};
    fitresultarray{i,3} = gof;
    
    % Extract rates and amplitudes from fittings and uncertainties
    % Extract rates
    fall_smallexp(i) = fitresult{1}.b;
    fall_mediumexp(i) = fitresult{1}.d;
    fall_largeexp(i) = fitresult{1}.g;
    rise_smallexp(i) = fitresult{2}.b;
    rise_largeexp(i) = fitresult{2}.d;
    
    %Extract amplitudes
    fall_smallamp(i) = fitresult{1}.a;
    fall_mediumamp(i) = fitresult{1}.c;
    fall_largeamp(i) = fitresult{1}.f;
    rise_smallamp(i) = fitresult{2}.a;
    rise_largeamp(i) = fitresult{2}.c;
    
    %Extract confidence bounds on all parameters
    confbnd_fall = confint(fitresult{1});
    confbnd_rise = confint(fitresult{2});
    
    %Extract uncertainties in rates
    fall_smallexp_lowbnd(i) = confbnd_fall(1,2);
    fall_mediumexp_lowbnd(i) = confbnd_fall(1,4);
    fall_largeexp_lowbnd(i) = confbnd_fall(1,6);
    rise_smallexp_lowbnd(i) = confbnd_rise(1,2);
    rise_largeexp_lowbnd(i) = confbnd_rise(1,4);
    
    fall_smallexp_error(i) = fall_smallexp(i) - fall_smallexp_lowbnd(i);
    fall_mediumexp_error(i) = fall_mediumexp(i) - fall_mediumexp_lowbnd(i);
    fall_largeexp_error(i) = fall_largeexp(i) - fall_largeexp_lowbnd(i);
    rise_smallexp_error(i) = rise_smallexp(i) - rise_smallexp_lowbnd(i);
    rise_largeexp_error(i) = rise_largeexp(i) - rise_largeexp_lowbnd(i);
    
    %Extract uncertainties in amplitudes
    fall_smallamp_lowbnd(i) = confbnd_fall(1,1);
    fall_mediumamp_lowbnd(i) = confbnd_fall(1,3);
    fall_largeamp_lowbnd(i) = confbnd_fall(1,5);
    rise_smallamp_lowbnd(i) = confbnd_rise(1,1);
    rise_largeamp_lowbnd(i) = confbnd_rise(1,3);
    
    fall_smallamp_error(i) = fall_smallamp(i) - fall_smallamp_lowbnd(i);
    fall_mediumamp_error(i) = fall_mediumamp(i) - fall_mediumamp_lowbnd(i);
    fall_largeamp_error(i) = fall_largeamp(i) - fall_largeamp_lowbnd(i);
    rise_smallamp_error(i) = rise_smallamp(i) - rise_smallamp_lowbnd(i);
    rise_largeamp_error(i) = rise_largeamp(i) - rise_largeamp_lowbnd(i);
end

% Save everything
savefilename = [filename_prefix,'.mat'];
save([savepath,savefilename], 'fitresultarray', 'fall_smallexp', 'fall_mediumexp', 'fall_largeexp', 'rise_smallexp', 'rise_largeexp', 'fall_smallexp_error', 'fall_mediumexp_error', 'fall_largeexp_error', 'rise_smallexp_error', 'rise_largeexp_error', 'fall_smallamp', 'fall_mediumamp', 'fall_largeamp', 'rise_smallamp', 'rise_largeamp', 'fall_smallamp_error', 'fall_mediumamp_error', 'fall_largeamp_error', 'rise_smallamp_error', 'rise_largeamp_error');







