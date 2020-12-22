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

 
    % Define the Lorentzian function.
    Lorentz = @(A,gamma,a0,B,x)...
        (A/pi)*(0.5*gamma)*(((x-a0).^2+0.25*(gamma^2)).^-1)+B;

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
xdata = [];
ydata = [];
hold on;

 for k = 1:N 
     
  Data1 = (data{k,1}.x(:));
  Data2 = (data{k,1}.y(:));
  xdata = Data1;
  ydata = Data2;
  
   % Here are the automated guesses for the 4 parameters, the first two are
        % straightforward. For the guess for gamma the FWHM is approximated based
        % on the data.
    
        % Added (6/14/18): Changed the guess for A1 by simply solving the equation
        % analytically like a normal person.
        a1 = xdata(ydata == max(ydata));
        B1 = min(ydata);
        
        halfmax = 0.5*(max(ydata)+min(ydata));
        index1 = find(ydata >= halfmax, 1, 'first');
        index2 = find(ydata >= halfmax, 1, 'last');
        gamma1 = xdata(index2)-xdata(index1);
        
        A1 = 0.5*pi*gamma1*(max(ydata)-B1);
        
        % Fit with a single Lorentzian.
        [f,gof] = fit(xdata,ydata,Lorentz,'StartPoint',[A1, gamma1, a1, B1]);
  
 end
 
 hold off;
  
  


