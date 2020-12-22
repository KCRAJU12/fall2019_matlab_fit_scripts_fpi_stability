

%% The script allows the user to select multiple .trc files once at a time.
%% These .trc files are read out by the function named as the ReadLecroyBinaryWaveform.m
%% function. The for loop is used to read all the data turn by turn.

[FileName,PathName] = uigetfile('*.trc', 'Open text file','MultiSelect','on');  % Helps to select required number of files
 data = cell(length(FileName),1) ;         % Saves the loaded data in the folder named "data" in the cell array form
 
 %% The for loop to read all the .trc file

 for i = 1:length(FileName) 
     file = ReadLeCroyBinaryWaveform(fullfile(PathName,FileName{i}));
     data{i} = file ;  % The curly brace for the for loop index is used if the data are saved in the cell array
 end

 

N = length(FileName);   % is length of the file uploaded by uigetfile function


%% Initialize of the Array of double data type to store Lorentzian fit coefficients

C=zeros(4);  % The double named C to store the fit coefficients from Lorentz fit in the double form

%% define the fitting function with fittype

%f=fittype(@(A,x0,gamma, x) (A./(pi*gamma))*(1+((x-x0).^2)./gamma).^-1);
 %Lorentz=fittype(@(A,x0,gamma,B,x)(A/pi)*(0.5*gamma)*(((x-x0).^2+0.50*(gamma^2)).^-1)+B); % Matlab automatically treats x as independent variable
% This Lorentzian is the area normalized Lorentzian. I have to use the
% amplitude normalized Lorentzian

Lorentz = fittype(@(A,x0,gamma,B,x) (A * (0.5 * gamma) * (((x-x0).^2+(0.5 * gamma).^2)).^-1)+B);

hold on;    % Helps to prevent the override of the file in each iteration of the for loop 

 for k = 1:N
   
  xdata = (data{k,1}.x(:));    % x column of the oscilloscope data; i.e.;
  %the time axis for the .trc file
  
  %xdata = data{k,1}(:,1);     % x column of the Oscilloscope data for the .txt file
  
 ydata = (data{k,1}.y(:)); % y column of the oscilloscope data; i.e.; the
 % transmission axis for the .trc file
 
% ydata = data{k,1}(:,2);     % x column of the Oscilloscope data for the .txt file 
  
 % Here are the automated guesses for the 4 parameters, the first two are
 % straightforward. For the guess for gamma the FWHM is approximated based
 % on the data.
 
        %a1 = mean(xdata(ydata == max(ydata))); % Note on a1 : If ydata is not unique, find(ydata==max(ydata)) can find more than one element.
        %a1 = (xdata(ydata == max(ydata)));
        [yvalue,loc] = max(ydata);
        a1 = xdata(loc);
        B1 = min(ydata);
        
        halfmax = 0.5*(max(ydata)+min(ydata));
        index1 = find(ydata >= halfmax, 1, 'first');
        index2 = find(ydata >= halfmax, 1, 'last');
        gamma1 = xdata(index2)-xdata(index1);
        
        %A1 = 0.5*pi*gamma1*(max(ydata)-B1);
        A1 = (max(ydata)-B1);
   
%% assign initial guessed parameters

% [A, x0, gamma,B] they are in the order of the appearance
% in the above fit function definition
% pin=[7.276 0.02309 0.0007 0.005]; % Initial guesses to fit the Gaussian function
% pin=[7.595848895609379 0.02461 0.00071 -0.275836229324341]; % Initial guesses to fit the Gaussian function

%% Finally, we are ready to fit our data

%[f,gof] = fit(xdata,ydata,Lorentz,'StartPoint',[A1, gamma1, a1, B1]);
fitobject = fit(xdata,ydata, Lorentz, 'StartPoint', [A1, a1, gamma1, B1]);     % is the fit object
%fitobject = fit(xdata,ydata, Lorentz, 'StartPoint', [6,0.02293,0.0013,0.06581]);     % is the fit object

C(k,:) = (coeffvalues(fitobject));                                             % is the matrix for the coefficients from the fit

%% Let's see how well our fit follows the data
figure;

plot(fitobject, xdata,ydata,'.g', 'fit')
%set(gca,'FontSize',24);                  % adjusting font size
%xlabel('Time(seconds)');                 % labels the x axis
%ylabel('Trans(volts)');                  %labels the y axis

 end
 
hold off;  % stops the loop after the process is complete

