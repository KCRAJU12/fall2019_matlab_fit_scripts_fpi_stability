 %% Fitting with Lorentzian(s) and extracting parameters.
% This program is used for loading and analyzing a bunch of RPLE spectra.
% Specfically comparing reference and AC scans to determine properties of
% the AC Stark shift and Dynamic Nuclear Polarization (DNP).

% ******* ver6.0 (12/18/18) **********
% UPDATES:
%
% The program now asks the user for a usertolerance value, so the user
% will not have to search through the program to change it.
%
% I got rid of the for loop for selecting a best value for
% 'MinPeakProminence' from previous versions, as it never really worked.
%
% Reference mean now prints to command line once it is calculated.
%
% All plotting arrays are now saved in the plotarrays structure. This now
% includes arrays for all 5 possible fitting parameters: w0s, linewidths,
% heights, areas, and background.
%
%
% Minor Update (8/7/19): the values for peak height had not been background
% subtracted in the initial rewrite, so I have now corrected this. There
% were also a few small errors in the manual fitting loop that I have now
% also fixed. The weights were also done improperely, so I have fixed this.

% PROGRAM OVERVIEW:  (goes roughly section by section)
% The program asks the user via uigetdir() to specify the
% path to the folder DIRECTLY ABOVE the folders with individual scan data.
% You can change the default folder this search starts in on line 82.
%
% The program loads all of the spectra into the "data" structure. It does
% this by searching all of the subfolders for a data file with the string
% 'spectrum' in the name.
%
% The program identifies each scan as either reference or AC by reading the
% name of the folders.
%
% The program determines if it has been run before for this set of data
% by searching the directory for the RPLE Analysis Data.m data file.
%
% The user is asked to enter a value for the tolerance the program will use
% to determine if fits are good enough. The user is also asked if this is
% important data that they would like to fit manually. If the user says no
% the program runs as normal.
%
% Each spectrum is fit with a single Lorentzian, then the gof.rsquared
% value is evaluated for AC scans to determine if fitting with a sum of two
% Lorentzians is warrented. This is dependent on detection polarization; X
% and Y are fit with a double, R and L are fit witht a single.
%
% If the user has indicated that this is important data, then the program
% asks the user to evaluate each fit (excluding references) by eye. If the
% user says the fit is good then the program continues as normal. If the
% user says the fit is bad then an interface to enter new guesses appears
% (the previous guesses are displayed as well). A plot of the new fit with
% the new guess is shown. This loop repeats until the user says the fit is
% good, then the program continues as normal.
%
% A structure called plotarrays is created for easy of plotting. Each
% array can be specified in the following way:
% plotarrays.fitparameter.(ref or AC)
%
% The user is asked if they would like to manually input some data to plot
% w0s against. For example power or polarization of AC laser. If the
% incorrect amount of data is entered, the user is given unlimited tries to
% enter the correct amount, feedback is provided. If the program has been
% run before the defualt input for userdata will be the array previously
% entered.
%
% Plots of (w0s,linewidths,heights,areas) vs. userdata are made with and
% without shaded regions depicting the FWHM of the peaks for the w0s plot.
%
% The figure and relevant variables are saved.
% ***********************************
 
 
    % Create an array called 'fitvalues' for all of the fit parameters.
    fitvalues = struct('filename',cellstr((chars(1:length(directory)-2))'),...
            'w0s',cellstr((chars(1:length(directory)-2))'),...
            'linewidths',cellstr((chars(1:length(directory)-2))'),...
            'heights',cellstr((chars(1:length(directory)-2))'),...
            'areas',cellstr((chars(1:length(directory)-2))'),...
            'B',cellstr((chars(1:length(directory)-2))'));
    
    % Define the Lorentzian function.
    Lorentz = @(A,gamma,a0,B,x)...
        (A/pi)*(0.5*gamma)*(((x-a0).^2+0.25*(gamma^2)).^-1)+B;
    
     % The for loop to fit all of the spectra.
    for i = 1:N
        fitvalues(i).filename = data(i).filename;
        
        xdata = data(i).x;
        ydata = data(i).y;
        yerr = data(i).sy;
        
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
        [f,gof] = fit(xdata,ydata,Lorentz,'StartPoint',[A1, gamma1, a1, B1],'Weights',max(yerr)./yerr);
        
         % Now output all of the fit parameters into 'fitvalues' structure.
        % If a single Lorentzian is used saves like [value, confint1, confint2]
        % If a double Lorentzian is used saves like
        % [value1, confint11, confint12, value2, confint21, confint22]
        values = coeffvalues(f);
        intervals = confint(f);
        
        if scantype(i) == 0 || scantype(i) == 1
            % Working out height from fitting parameters.
            h = (2*values(1))/(pi*values(2)) + values(4);
            hint1 = (2*intervals(1,1))/(pi*values(2)) + values(4);
            hint2 = (2*intervals(2,1))/(pi*values(2)) + values(4);
            
            % Place fit parameters in the fitvalues structure.
            fitvalues(i).heights = [ h hint1 hint2 ];
            fitvalues(i).linewidths = [ values(2) intervals(1,2) intervals(2,2) ];
            fitvalues(i).w0s = [ values(3) intervals(1,3) intervals(2,3) ];
            fitvalues(i).B = [ values(4) intervals(1,4) intervals(2,4) ];
        elseif scantype(i) == 2
            % Working out height from fitting parameters.
            h1 = (2*values(1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
            h1int1 = (2*intervals(1,1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
            h1int2 = (2*intervals(2,1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
            
            h2 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*values(2))/(pi*values(4)) + values(7);
            h2int1 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*intervals(1,2))/(pi*values(4)) + values(7);
            h2int2 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*intervals(2,2))/(pi*values(4)) + values(7);
            
            
            fitvalues(i).heights = [ h1 h1int1 h1int2...
                h2 h2int1 h2int2 ];
            
            fitvalues(i).areas = [ values(1) intervals(1,1) intervals(2,1)...
                values(2) intervals(1,2) intervals(2,2) ];
            
            fitvalues(i).linewidths = [ values(3) intervals(1,3) intervals(2,3)...
                values(4) intervals(1,4) intervals(2,4) ];
            
            fitvalues(i).w0s = [ values(5) intervals(1,5) intervals(2,5)...
                values(6) intervals(1,6) intervals(2,6) ];
            
            fitvalues(i).B = [ values(7) intervals(1,7) intervals(2,7) ];
        end
    end
    % Delete parts of 'fitvalues' that were unused.
    for i = length(fitvalues):-1:1
        if ischar(fitvalues(i).w0s) == 1
            fitvalues(i)=[];
        end
    end
    
    % Check to confirm that data were loaded and fit appropriately.
    if N ~= length(fitvalues)
        cprintf('err', '\nERROR: The length of the data array does not match the length of the fitvalues array.\n');
        beep; return
    end 
    
    
    
    