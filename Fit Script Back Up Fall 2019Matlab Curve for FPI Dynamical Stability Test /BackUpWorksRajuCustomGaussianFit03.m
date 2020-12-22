%% The script allows the user to select multiple .trc files once at a time.
%% These .trc files are read out by the function named as the ReadLecroyBinaryWaveform.m
%% function. The for loop is used to read all the data turn by turn.

[FileName,PathName] = uigetfile('*.trc', 'Open text file','MultiSelect','on');  % Helps to select required number of files
 data = cell(length(FileName),1) ;         % Saves the loaded data in the folder named "data" in the cell array form
 
 %% The for loop to read all the .trc file

 for i = 1:length(FileName) 
     file = ReadLeCroyBinaryWaveform(fullfile(PathName,FileName{i}));
     data{i} = file ;
 end
 
 %% Initialization of the variable for the fitting purpose
 
 x = [];  % Initialization of the independent variable
 y = [];  % Initialization of the dependent variable
 
tic;     % Is the starting time of the loop

N = length(FileName);   % is length of the file uploaded by uigetfile function


%% Initialize of the Array of double data type to store Lorentzian fit coefficients

C=zeros(4);  % The double named C to store the fit coefficients from Lorentz fit in the double form

hold on;    % Helps to prevent the override of the file in each iteration of the for loop 
 for k = 1:N
   
  Data1 = (data{k,1}.x(:));
  Data2 = (data{k,1}.y(:));
  x = Data1;
  y = Data2;
   
%% define the fitting function with fittype

%f=fittype(@(A,x0,gamma, x) (A./(pi*gamma))*(1+((x-x0).^2)./gamma).^-1);
 f=fittype(@(A,x0,gamma,B,x)(A/pi)*(0.5*gamma)*(((x-x0).^2+0.50*(gamma^2)).^-1)+B); % Matlab automatically treats x as independent variable


%% assign initial guessed parameters

% [A, x0, gamma] they are in the order of the appearance
% in the above fit function definition
%pin=[7.276 0.02309 0.0007 0.005]; % Initial guesses to fit the Gaussian function
pin=[7.403856575489044 0.02461 0.00072 -0.275836229324341]; % Initial guesses to fit the Gaussian function

%% Finally, we are ready to fit our data

[fitobject] = fit (x,y, f, 'StartPoint', pin);
C(k,:)=(coeffvalues(fitobject));

%% Let's see how well our fit follows the data

plot(fitobject, x,y,'-k', 'fit')
set(gca,'FontSize',24);                  % adjusting font size
xlabel('Time(seconds)');                 % labels the x axis
ylabel('Trans(volts)');                  %labels the y axis


 end
 
hold off;  % stops the loop after the process is completed

toc; % determines the elapsed time

