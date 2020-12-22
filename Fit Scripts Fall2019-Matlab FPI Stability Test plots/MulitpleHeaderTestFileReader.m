
%% The Oscilloscope saves the data in the text format with multiple number of heading 
%% The script written below reads these headers in the text files and save the data in the double format with xdata
%% named as "xdata" and the ydata named as "ydata"

[FileName,PathName] = uigetfile('*.txt', 'Open text file','MultiSelect','on');  % Helps to select required number of files
 data = cell(length(FileName),1) ;         % Saves the loaded data in the folder named "data" in the cell array form
 
 %% The for loop to read all the .trc file

 for i = 1:length(FileName) 
     file = importdata(fullfile(PathName,FileName{i}));
     data{i} = file ;  % The curly brace for the for loop index is used if the data are saved in the cell array
     %xdata = str2double(data{i,1}.textdata(6:end,1));
     %ydata = str2double(data{i,1}.textdata(6:end,2));
 end

 % The following script read a single text file with the multiple number of header
%% lines Keep it for the future purpose.Do not delete the commented section. If you need, we can 
%% Uncomment it and use the command importdata to read the multiple header text file
%A = importdata('Z3Trace00002.txt');
%x = A.textdata(6:end,1);
%y = A.textdata(6:end,2);
%new = str2double(x); 
 
tic;                    % Is the starting time of the loop

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
   
   xdata = str2double(data{k,1}.textdata(6:end,1));% x column of the oscilloscope data; i.e.; the time axis
   ydata = str2double(data{k,1}.textdata(6:end,2)); % y column of the oscilloscope data; i.e.; the transmission axis
  
 % Here are the automated guesses for the 4 parameters, the first two are
 % straightforward. For the guess for gamma the FWHM is approximated based
 % on the data.
 
        %a1 = mean(xdata(ydata == max(ydata))); % Note on a1 : If ydata is not unique, find(ydata==max(ydata)) can find more than one element.
        %a1 = (xdata(ydata == max(ydata)));
        [yvalue,loc] = max(ydata); % locates the maximum ydata value
        a1 = xdata(loc); % Corresponding peak center related to maximum yadata along the xaxis
        B1 = min(ydata); % Minimum ydata value
        
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
%fitobject = fit(xdata,ydata, Lorentz, 'StartPoint', [A1, gamma1, a1, B1]);     % is the fit object
fitobject = fit(xdata,ydata, Lorentz, 'StartPoint', [A1,a1,gamma1,B1]);     % is the fit object

C(k,:) = (coeffvalues(fitobject));                                             % is the matrix for the coefficients from the fit

%% Let's see how well our fit follows the data
figure;

plot(fitobject, xdata,ydata,'.g', 'fit')
set(gca,'FontSize',24);                  % adjusting font size
xlabel('Time(seconds)');                 % labels the x axis
ylabel('Trans(volts)');                  %labels the y axis

 end
 
hold off;  % stops the loop after the process is completed

toc; % determines the elapsed time
