%% Upload some data which will be used as the dependent and the independent 
%variable in the Matlab.

[FileName,PathName] = uigetfile('*.trc', 'Open text file','MultiSelect','on');
 data = cell(length(FileName),1) ;
 for i = 1:length(FileName) 
     file = ReadLeCroyBinaryWaveform(fullfile(PathName,FileName{i}));
     data{i} = file ;
 end
 
 x = []; % Independent variable
 y = []; % Dependent variable
tic;

N = length(FileName); % is length of the file uploaded by uigetfile function

%% Initialize Arrays for Extracted Lorentzian fit coefficients
Amp_Coeff = zeros(N,1);      % is amplitude in the Lorentz fit function
PeakCtr_Coeff = zeros(N,1);  % is the peak center of the fitted Lorentzian
FWHM_Coeff = zeros(N,1);     % is the Full widtha at Half Maximum of fitted Lorentzian
Extra_Coeff = zeros(N,1);    % is the adjustment coefficient in the Lorentzian

fitresultarray = cell(N,3); % Contains all information for the Lorentzian fits


% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 2, 1 );
gof = struct( 'sse', cell( 2, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
hold on;
 for k = 1:N
   
  Data1 = (data{k,1}.x(:));
  Data2 = (data{k,1}.y(:));
  x = Data1;
  y = Data2;
   
%% define the fitting function with fittype
% Matlab automatically treats x as independent variable
%f=fittype(@(A,x0,gamma, x) (A./(pi*gamma))*(1+((x-x0).^2)./gamma).^-1);
% f=fittype(@(A,x0,gamma,B,x)(A/pi)*(0.5*gamma)*(((x-x0).^2+0.25*(gamma^2)).^-1)+B);
 f=fittype(@(A,x0,gamma,B,x)(A/pi)*(0.5*gamma)*(((x-x0).^2+0.50*(gamma^2)).^-1)+B);


%% assign initial guessed parameters
% [A, x0, gamma] they are in the order of the appearance
% in the above fit function definition
%pin=[7.276 0.02309 0.0007 0.005]; % Initial guesses to fit the Gaussian function
pin=[7.403856575489044 0.02461 0.00072 -0.275836229324341]; % Initial guesses to fit the Gaussian function

%% Finally, we are ready to fit our data

[fitobject] = fit (x,y, f, 'StartPoint', pin);


%% Let's see how well our fit follows the data
plot(fitobject, x,y,'-k', 'fit')
set(gca,'FontSize',24); % adjusting font size
xlabel('Time(seconds)');
ylabel('Trans(volts)');

%% Establish arrays to contain all the fit parameters
fitresultarray{k,1} = fitresult{1};

%% Extract the coefficients from the Lorentz fit

Amp_Coeff(k) = fitresult{1}.A;          % is amplitude in the Lorentz fit function

PeakCtr_Coeff(k) = fitresult{1}.x0;     % is the peak center of the fitted Lorentzian

FWHM_Coeff(k) = fitresult{1}.gamma;     % is the Full widtha at Half Maximum of fitted Lorentzian

Extra_Coeff(k) = fitresult{1}.B;        % is the adjustment coefficient in the Lorentzian

 end
hold off;

toc;

