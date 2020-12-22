% This script loads the multiple text files, reads them, does the
% gaussian fit to the data and saves the extracted coefficients from the
% fit in the cell array format
%[filename, pathname, filterindex] = uigetfile({'*.mat',  'Matlab data files (*.mat)'; '*.*', 'All Files (*.*)'}, 'Select a data set', 'MultiSelect', 'on');  % loads the files
[filename, pathname, filterindex] = uigetfile('*.mat','MultiSelect','on');
N = length(filename);                                                     % Number of files loaded


% Now I will initialize the Arrays for the Extracted exponentials from the
% fit

%Initialize arrays for extracted exponentials from single gauss1 fit

amp_coeff1 = zeros(N,1);    %Amplitude coefficeints array for height of the bell curve in gauss1 function
centr_coeff1 = zeros(N,1);  %Center coeeficients array for the center of the gaussian peaks
stdv_coeff1 = zeros(N,1);   %Standard deviation coefficients array for the gaussian peaks

%amp_coeff2 = zeros(N,1);    %Amplitude coefficeints array for height of the bell curve in gauss1 function
%centr_coeff2 = zeros(N,1);  %Center coeeficients array for the center of the gaussian peaks
%stdv_coeff2 = zeros(N,1);   %Standard deviation coefficients array for the gaussian peaks

fitresultarray = cell(N,3);  % This contains all informations for fit

for i = 1:N
load([pathname,filename{i}]);
    
    % Fit current data
    
    
[fitresult, gof] = createFit(X, Y);      % The calling of the
    

    
%initialguess = [0.1427  0.0648 0.0200];            % Samples data points for the fit

                                              
%Establish arrays to contain all the parameters
fitresultarray{i,1} = fitresult{1};
%fitresultarray{i,2} = fitresult{2};
fitresultarray{i,3} = gof;

amp_coeff1 =   fitresult{1}.a1;    % Extracted array of amplitude coefficients from fit
centr_coeff1 = fitresult{1}.a2;    % Extracted array of Gaussian's peak coefficients from fit
stdv_coeff1 =  fitresult{1}.a3;    % Extracted array of standard deviation coefficients from fit

%amp_coeff2 =   fitresult{1}.a2;    % Extracted array of amplitude coefficients from fit
%centr_coeff2 = fitresult{1}.b2;    % Extracted array of Gaussian's peak coefficients from fit
%stdv_coeff2 =  fitresult{1}.c2;    % Extracted array of standard deviation coefficients from fit


end

X = A(:,1);
Y = A(:,2);










