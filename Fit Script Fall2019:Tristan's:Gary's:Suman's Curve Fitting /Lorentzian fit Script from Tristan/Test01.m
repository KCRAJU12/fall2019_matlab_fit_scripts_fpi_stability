function RPLE_Analysis_ver6_0
%% RPLE_Analysis_ver6_0.m
% By Tristan Wilkinson

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
%
%
%
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
% The program determines if the it has been run before for this set of data
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

%% Loading in the Data
clearvars;
close all;
clc;

% Ask the user to select the folder of interest. Start point can be
% changed by specifying the string inside of uigetdir().
path = uigetdir('\\ecas.wvu.edu\squol\AC Stark Effect\Charge Control Sample\AC Stark Paper Data\z_Refitting new PL mask');

if path == 0 %User pressed cancel
    cprintf('err', '\nERROR: Please select the folder DIRECTLY ABOVE the folders with data in them when prompted.\n');
    beep; return
end

% Extract the name of the folder the user has selected.
split = strsplit(path, '\');
folder = split(end);

fprintf(1, ['\nRPLE_Analysis: ', folder{1}, '\n']);

% Load in the structure for the folder the program is in.
directory = dir(path);

% Next define a 'data' structure to store the filename,x,y,sy in when
% the spectrum is loaded in. (Preallocation makes code run faster)
chars = char(1:length(directory));

data = struct('filename',cellstr((chars(1:length(directory)-2))'),...
    'x',cellstr((chars(1:length(directory)-2))'),...
    'y',cellstr((chars(1:length(directory)-2))'),...
    'sy',cellstr((chars(1:length(directory)-2))'));

% Define the value of speed of light for conversions to frequency.
    c = 299792458;
    
    % For loop to load in the spectrum and assign them to the 'data'
    % structure sequentially. Because there may be other files in this folder
    % there will be extra spaces in the structure.
    for i = 3:length(directory)
        if directory(i).isdir == 1
            subfolder = dir([path '\' directory(i).name]);
            spectrumfound = 0; %Local feedback
            for j = 3:length(subfolder) %Search the subfolder for any files that have 'spectrum' in them
                split1 = strsplit(subfolder(j).name, ' ');
                isspectrumA = ismember(split1, 'spectrum');
                isspectrumB = ismember(split1, 'Spectrum');
                if any( [isspectrumA isspectrumB] )
                    k = j;
                    spectrumfound = 1;
                end
            end
            if spectrumfound == 0 %If no spectrum file is found
                cprintf('err', ['\nERROR: There is no spectrum file in ' directory(i).name ' .\n']);
                beep; return
            end
            load([ path '\' directory(i).name '\' subfolder(k).name ], 'x', 'y', 'sy');
            data(i-2).filename = directory(i).name; %Loads in the filename
            data(i-2).x = c./x; %Convert wavelength to frequency
            data(i-2).y = y;
            data(i-2).sy = sy;
        else
            continue,
        end
    end
    
    % This loop goes through in reverse order to delete the parts of the 'data'
    % structure that were left unused.
    for i = length(data):-1:1
        if ischar(data(i).x) == 1
            data(i)=[];
        end
    end
    
    % Define number of total scans as N.
    N = length(data);
    
    fprintf(1, ['\n' num2str(N) ' RPLE spectra loaded in:\n']);
    
    %% Determining if scans are reference or AC.
    % This section parses the filename strings to look for a string of
    % 'reference' or 'AC'.
    
    % Array for scan type, 0 means reference and 1 means AC.
    scantype = zeros(1, N);
    
    for i = 1:N
       split2 = strsplit(data(i).filename, ' ');
       refscanA = ismember(split2,'reference');
       refscanB = ismember(split2,'Reference');
       ACscanA = ismember(split2,'AC');
       ACscanB = ismember(split2,'ac');
       if any( [refscanA refscanB] )
           scantype(i) = 0;
       elseif any( [ACscanA ACscanB] )
           scantype(i) = 1;
       else
           cprintf('err', ['\nERROR: Folder ' data(i).filename ' is not specified as reference or AC.\n']);
           beep; return
       end
    end
    
    % Determine the number of reference and AC scans.
    nref = ismember(scantype,0);
    nAC = ismember(scantype,1);
    
    Nref = sum(nref);
    NAC = sum(nAC);
    
    fprintf(1, ['\t\t\t\t\t\t\t', num2str(Nref), ' reference\n']);
    fprintf(1, ['\t\t\t\t\t\t\t', num2str(NAC), ' AC\n']);
    
    %% Determine if program ran previously.
    warning('off','all') %Turn warnings off so matlab does not spit one if the variables below are not found
    
    ranbefore = 0;
    userenteredbefore = 0;
    for i = 3:length(directory)
        split2a = strsplit(directory(i).name, ' ');
        if length(split2a) > 1 %If there are no spaces split2a{end-1}=0 and there is an error, so avoid this
            if strcmp([split2a{end-1} ' ' split2a{end}], 'RPLE data.mat')
                ranbefore = 1;
                load([ path '\' directory(i).name ], 'usertolerance', 'userdata', 'userdatalabel');
                usertolenteredbefore = exist('usertolerance'); %Incase older version of program has been ran before that did not save usertolerance
                userenteredbefore = exist('userdata'); %Incase program has been ran before but user chose not to enter data
            else
                continue;
            end
        end
    end
    
    %% Fitting with Lorentzian(s) and extracting parameters.
    
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
    
    % Define the sum of two Lorentzians as a function.
    Lorentz2 = @(A,A2,gamma,gamma2,a0,a02,B,x)...
        (A/pi)*(0.5*gamma)*(((x-a0).^2+0.25*(gamma^2)).^-1)+...
        (A2/pi)*(0.5*gamma2)*(((x-a02).^2+0.25*(gamma2^2)).^-1)+B;
    
    % Define a global feedback variable to know if any scans are fit with
    % the sum of two Lorentzians. And if any scans could not be fit well.
    feedbackglobal = 0;
    
    % Set a value for tolerance used to decide how many peaks to fit with.
    % If the program has been ran before provide option to choose the same.
    % Two different dialog boxes may appear based on if ranbefore.
    if ranbefore == 0 || (ranbefore == 1 && usertolenteredbefore == 0)
        answertol = questdlg(sprintf([ 'What value for user tolerance would you like to use?'...
            '\n\n0.92 is default'...
            '\n0.6 to fit with one peak (Linear AC or only R/L detection)' ]),...
            'User tolerance value entry',...
            '0.92', '0.60', 'Custom', '0.92');
    elseif ranbefore == 1 && usertolenteredbefore == 1
        answertol = questdlg(sprintf([ 'What value for user tolerance would you like to use?'...
            '\n\n0.92 is default'...
            '\n0.6 to fit with one peak (Linear AC or only R/L detection)'...
            '\n' num2str(usertolerance) ' is previously entered value' ]),...
            'User tolerance value entry',...
            '0.92', num2str(usertolerance), 'Custom', num2str(usertolerance));
    end
    
    % Assign the value the user selected.
    switch answertol
        case '0.92'
            usertolerance = 0.92;
        case '0.60'
            usertolerance = 0.6;
        case '' %User closed the dialog box
            cprintf('err', '\nERROR: Please enter a user tolerance value.\n');
            beep; return
        case 'Custom' %Let the user enter their desired value manually
            answercustomtol = inputdlg('User tolerance value:',...
                'Custom user tolerance value entry',...
                [1 50]);
            usertolerance = str2num(answercustomtol{1});
    end
    
    fprintf(1, '\nFitting progress:\n');
    
    % Ask the user whether this is important data that needs each fit to be
    % checked manually by the user. (Basically only necessary if this is
    % final data for a paper)
    answer0 = questdlg('Would you like to fit this data manually?',...
        'Manual fitting');
    
    switch answer0
        case 'No'
            manualfit = 0;
        case 'Cancel'
            cprintf('err', '\nERROR: Please select either yes or no when asked if you would like to fit scans manually.\n');
            beep; return
        case '' %User closed the dialog box
            cprintf('err', '\nERROR: Please select either yes or no when asked if you would like to fit scans manually.\n');
            beep; return
        case 'Yes' %This is important data the user would like to check each fit manually.
            manualfit = 1;
    end
    
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
        a1 = xdata(find(ydata == max(ydata)));
        B1 = min(ydata);
        
        halfmax = 0.5*(max(ydata)+min(ydata));
        index1 = find(ydata >= halfmax, 1, 'first');
        index2 = find(ydata >= halfmax, 1, 'last');
        gamma1 = xdata(index2)-xdata(index1);
        
        A1 = (max(ydata)-B1);
        
        % Fit with a single Lorentzian.
        [f,gof] = fit(xdata,ydata,Lorentz,'StartPoint',[A1, gamma1, a1, B1],'Weights',max(yerr)./yerr);
        
        % Evaluate the fits rsquared value to see if fitting with a sum of
        % two Lorentzians is warranted.
        feedbacklocal = 0;
        badfit = 0;
        
        % Determine whether to fit with two or one Lorentzians based on
        % detection polarization.
        split3 = strsplit(data(i).filename, ' ');
        X = ismember(split3,'X');
        Y = ismember(split3,'Y');
        R = ismember(split3,'R');
        L = ismember(split3,'L');
        if any(X) || any(Y) %Fit X and Y detection with 2 peaks
            tolerance = usertolerance;
        elseif any(R) || any(L) %Fit R and L detection with 1 peak
            tolerance = 0.6;
        else
            tolerance = usertolerance;
        end
        
        if scantype(i) == 1 && gof.rsquare <= (tolerance + 0.077)
            feedbacklocal = 1; %Local knowledge of if statement outcome
            feedbackglobal = 1; %Global knowledge of if statement outcome
            
            minPeak = 12000;
            [t,q] = findpeaks(ydata, 'NPeaks', 2, 'MinPeakProminence', minPeak);
            
            j = 0;
            while isempty(t)
                j = j+100;
               [t,q] = findpeaks(ydata, 'NPeaks', 2, 'MinPeakProminence', minPeak-j);
            end
            
            % Fit the data with a sum of two Lorentzians.
            if length(q) == 1 %findpeaks() finds only a single peak
                q1 = [ q-1 q-1 ];
                t1 = [ t t ];
                
                a2 = xdata(q1(1));
                a3 = xdata(q1(2)) + gamma1/4;
                gamma2 = gamma1/2;
                
                % These guesses come from plugging in the two peak points into
                % the sum of two Lorentzians, and then solving.
                A2 = (pi*gamma2*(4*(a2-a3).^2+gamma2.^2)*(-4*(B1-t1(1))*(a2-a3).^2+(t1(1)-t1(2))*gamma2.^2))/(...
                    8*(a2-a3).^2*(4*(a2-a3).^2+gamma2.^2)+8*(a2-a3).^2*gamma2.^2);
                A3 = -(pi*(4*(B1-t1(2))*(a2-a3).^2+(t1(1)-t1(2))*gamma2.^2)*gamma2*(4*(a2-a3).^2+gamma2.^2))/(...
                    8*(a2-a3).^2*(4*(a2-a3).^2+gamma2^2)+8*(a2-a3).^2*gamma2.^2);
                
                %Fit with a sum of two Lorentzians.
                [f,gof] = fit(xdata,ydata,Lorentz2,'StartPoint',[A2, A3, gamma2, gamma2, a2, a3, B1],'Weights',max(yerr)./yerr);
                
            elseif length(q) == 2 %findpeaks() finds two peaks.
                q2 = [ q(1)-1 q(2)-1 ];

                a2 = xdata(q2(1));
                a3 = xdata(q2(2));
                gamma2 = gamma1/5;
                
                % These guesses come from pluggin in the two peak points into
                % the sum of two Lorentzians, and then solving.
                A2 = (pi*gamma2*(4*(a2-a3).^2+gamma2.^2)*(-4*(B1-t(1))*(a2-a3).^2+(t(1)-t(2))*gamma2.^2))/(...
                    8*(a2-a3).^2*(4*(a2-a3).^2+gamma2.^2)+8*(a2-a3).^2*gamma2.^2);
                A3 = -(pi*(4*(B1-t(2))*(a2-a3).^2+(t(1)-t(2))*gamma2.^2)*gamma2*(4*(a2-a3).^2+gamma2.^2))/(...
                    8*(a2-a3).^2*(4*(a2-a3).^2+gamma2^2)+8*(a2-a3).^2*gamma2.^2);
                
                %Fit with a sum of two Lorentzians.
                [f,gof] = fit(xdata,ydata,Lorentz2,'StartPoint',[A2, A3, gamma2, gamma2, a2, a3, B1],'Weights',max(yerr)./yerr);
            end
            
            % Check to see if the fit is good.
            if gof.rsquare <= tolerance
                badfit = 1;
                fprintf(1, ['\t\t\t\t' data(i).filename ' has not been fit well\n']);
            end
            
            % Set the scantype array to be: 0 for reference, 1 for AC with
            % 1 peak, 2 for AC with 2 peaks. This effectively replaces
            % numpksAC from previous versions of program.
            if feedbacklocal == 1
                scantype(i) = 2;
            end
        end
        
        % Plot the raw data with fits.
        figure('Name', data(i).filename,'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h0 = errorbar(xdata,ydata,yerr,'b.','Capsize',0.1);
        h1 = plot(f);
        
        % If the fit is bad, show the guessed Lorentzians.
        if badfit == 1 && length(q) == 1
            h2 = plot(xdata, Lorentz2(A2,A3,gamma2,gamma2,a2,a3,B1,xdata), 'k'); %Guess will appear black
            legend([h0,h1,h2],'data', 'fit', 'guess')
        elseif badfit == 1 && length(q) == 2
            h2 = plot(xdata, Lorentz2(A2,A3,gamma2,gamma2,a2,a3,B1,xdata), 'g'); %Guess will appear green
            legend([h0,h1,h2],'data', 'fit', 'guess')
        else
            legend([h0,h1],'data','fit')
        end
        
        xlabel('Frequency (GHz)');
        ylabel('Intensity (arb. units)');
        
        % If the user has specified that this is important data, allow them
        % to look at the fit and decide if it is good enough.
        if manualfit == 1 && ( scantype(i) == 1 || scantype(i) == 2 ) %Only worry about this for AC scans
            answerfit = questdlg('Is this fit okay?','Manual Fitting');
            
            switch answerfit
                case 'No' %The fit is bad and the user would like to enter guesses manually
                    manualbadfit = 1;
                    % Variables to tell how the user would like to refit
                    % the data.
                    singletosingle = 0;
                    singletodouble = 0;
                    doubletosingle = 0;
                    doubletodouble = 0;
                    
                    if scantype(i) == 1 %This scan fit with a single Lorentzian
                        % Ask if the user would like to fit with 2
                        % Lorentzians.
                        answersingle = questdlg('How many Lorentzians would you like to fit with?',...
                            'Manual Fitting',...
                            'One','Two','Two');
                        
                        switch answersingle
                            case 'One' %This seems unlikely, so I will table it unless it comes up
                                singletosingle = 1;
                            case ''
                                cprintf('err', '\nERROR: Please select either yes or no when asked how many Lorentzians to fit with.\n');
                                beep; return
                            case 'Two'
                                singletodouble = 1;
                                scantype(i) = 2;
                        end
                        
                    elseif scantype(i) == 2 %This scan fit with a double Lorentzian
                        % Ask if the user would like to fit with 2
                        % Lorentzians.
                        answerdouble = questdlg('How many Lorentzians would you like to fit with?',...
                            'Manual Fitting',...
                            'One','Two','Two');
                        
                        switch answerdouble
                            case 'One' %This seems unlikely, so I will table it unless it comes up
                                doubletosingle = 1;
                                scantype(i) = 1;
                            case ''
                                cprintf('err', '\nERROR: Please select either yes or no when asked how many Lorentzians to fit with.\n');
                                beep; return
                            case 'Two'
                                doubletodouble = 1;
                        end
                        
                    end
                case 'Cancel'
                    cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                    beep; return
                case '' %User closed the dialog box
                    cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                    beep; return
                case 'Yes' %The fit is good, so move on
                    manualbadfit = 0;
            end
            % Now we know how many Lorentzians the initial data was fit
            % with, and how many the user would like to manually fit with.
            % Next run a while loop until to refit until the user is
            % satisfied with the results.
            
            j = 1; %Index the while loop so we know if it is first guess or not
            while manualbadfit == 1 %Repeat until the user says the fit is good enough
                answerboxdim = 50;
                if singletosingle == 1 || doubletosingle == 1 %The user wants to fit with a single Lorentzian
                    if singletosingle == 1
                        prompt = {[ 'Guess for A (Previous was ' num2str(A1) '):' ],...
                            [ 'Guess for gamma (Previous was ' num2str(gamma1) '):' ],...
                            [ 'Guess for w0 (Previous was ' num2str(a1) '):' ]};
                        title = 'Manual Guess Entry';
                        dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                        definput = {num2str(A1),num2str(gamma1),num2str(a1)};
                        opts.Resize = 'on';
                        opts.WindowStyle = 'normal';

                        answersinglefit = inputdlg(prompt,title,dims,definput,opts);
                    
                    elseif doubletosingle == 1 && j == 1
                        prompt = {[ 'Guess for A (Previous were ' num2str(A2) ', ' num2str(A3) '):' ],...
                            [ 'Guess for gamma ' num2str(gamma2) ', ' num2str(gamma2) '):' ],...
                            [ 'Guess for w0 ' num2str(a2) ', ' num2str(a3) ')' ]};
                        title = 'Manual Guess Entry';
                        dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                        definput = {num2str(A2),num2str(gamma2),num2str(a2)};
                        opts.Resize = 'on';
                        opts.WindowStyle = 'normal';

                        answersinglefit = inputdlg(prompt,title,dims,definput,opts);
                        
                    elseif doubletosingle == 1 && j > 1
                        prompt = {[ 'Guess for A (Previous was ' num2str(A1) '):' ],...
                            [ 'Guess for gamma (Previous was ' num2str(gamma1) '):' ],...
                            [ 'Guess for w0 (Previous was ' num2str(a1) '):' ]};
                        title = 'Manual Guess Entry';
                        dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                        definput = {num2str(A1),num2str(gamma1),num2str(a1)};
                        opts.Resize = 'on';
                        opts.WindowStyle = 'normal';

                        answersinglefit = inputdlg(prompt,title,dims,definput,opts);
                    end
                    
                    % Refit single Lorentzian with user guesses.
                    a1 = str2num(answersinglefit{3});
                    B1 = min(ydata);
                    gamma1 = str2num(answersinglefit{2});
                    A1 = str2num(answersinglefit{1});
                    
                    [f,~] = fit(xdata,ydata,Lorentz,'StartPoint',[A1, gamma1, a1, B1],'Weights',max(yerr)./yerr);
                    
                    % Plot the fit again so the user can look at it.
                    clf %Clear the current figure window
                    hold on
                    h0 = errorbar(xdata,ydata,yerr,'b.','Capsize',0.1); %Plot in the same figure window, overwriting the old one
                    h1 = plot(f);
                    
                    % Plot the guess as well to aid the user in fitting
                    h2 = plot(xdata, Lorentz(A1,gamma1,a1,B1,xdata), 'g');
                    legend([h0,h1,h2],'data', 'fit', 'guess')
                    
                    % Ask again if the fit is good.
                    answerfit = questdlg('Is this fit okay?','Manual fitting');
                    
                    switch answerfit
                        case 'No'
                            manualbadfit = 1;
                        case 'Cancel'
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case ''
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case 'Yes'
                            scantype(i) = 1; %Keep track of number of peaks for analysis later on
                            
                            % Get rid of the guess on the plot
                            clf
                            h0 = errorbar(xdata,ydata,yerr,'b.','Capsize',0.1);
                            h1 = plot(f);
                            legend([h0,h1],'data', 'fit')
                            
                            manualbadfit = 0; %Exit the while loop
                    end

                elseif singletodouble == 1 || doubletodouble == 1 %The user wants to fit with double Lorentzian
                    if singletodouble == 1 && j == 1
                        prompt = {[ 'Guess for A1 (Previous was ' num2str(A1) '):' ],...
                            'Guess for A2 (Previous was NA):',...
                            [ 'Guess for gamma1 (Previous was ' num2str(gamma1) '):' ],...
                            'Guess for gamma2 (Previous was NA):',...
                            [ 'Guess for w01 (Previous was ' num2str(a1) '):' ],...
                            'Guess for w02 (Previous was NA):',...
                            [ 'Guess for B (Previous was ' num2str(B1) '):' ]};
                        title = 'Manual Guess Entry';
                        dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                        definput = {num2str(A1),num2str(A1),num2str(gamma1),num2str(gamma1),...
                            num2str(a1),num2str(a1),num2str(B1)};
                        opts.Resize = 'on';
                        opts.WindowStyle = 'normal';

                        answerdoublefit = inputdlg(prompt,title,dims,definput,opts);
                        
                    elseif singletodouble == 1 && j > 1
                        prompt = {[ 'Guess for A1 (Previous was ' num2str(A2) '):' ],...
                            [ 'Guess for A2 (Previous was ' num2str(A3) '):' ],...
                            [ 'Guess for gamma1 (Previous was ' num2str(gamma2) '):' ],...
                            [ 'Guess for gamma2 (Previous was ' num2str(gamma2) '):' ],...
                            [ 'Guess for w01 (Previous was ' num2str(a2) '):' ],...
                            [ 'Guess for w02 (Previous was ' num2str(a3) '):' ],...
                            [ 'Guess for B (Previous was ' num2str(B1) '):' ]};
                        title = 'Manual Guess Entry';
                        dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                        definput = {num2str(A2),num2str(A3),num2str(gamma2),num2str(gamma2),...
                            num2str(a2),num2str(a3),num2str(B1)};
                        opts.Resize = 'on';
                        opts.WindowStyle = 'normal';

                        answerdoublefit = inputdlg(prompt,title,dims,definput,opts);
                        
                    elseif doubletodouble == 1
                        if j ==  1
                           gamma3 = gamma2; 
                        end
                        
                        prompt = {[ 'Guess for A1 (Previous was ' num2str(A2) '):' ],...
                            [ 'Guess for A2 (Previous was ' num2str(A3) '):' ],...
                            [ 'Guess for gamma1 (Previous was ' num2str(gamma2) '):' ],...
                            [ 'Guess for gamma2 (Previous was ' num2str(gamma3) '):' ],...
                            [ 'Guess for w01 (Previous was ' num2str(a2) '):' ],...
                            [ 'Guess for w02 (Previous was ' num2str(a3) '):' ],...
                            [ 'Guess for B (Previous was ' num2str(B1) '):' ]};
                        title = 'Manual Guess Entry';
                        dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                        definput = {num2str(A2),num2str(A3),num2str(gamma2),num2str(gamma3),...
                            num2str(a2),num2str(a3),num2str(B1)};
                        opts.Resize = 'on';
                        opts.WindowStyle = 'normal';

                        answerdoublefit = inputdlg(prompt,title,dims,definput,opts);
                    end
                    
                    % Refit double Lorentzian with user guesses.
                    a2 = str2num(answerdoublefit{5});
                    a3 = str2num(answerdoublefit{6});
                    gamma2 = str2num(answerdoublefit{3});
                    gamma3 = str2num(answerdoublefit{4});
                    
                    A2 = str2num(answerdoublefit{1});
                    A3 = str2num(answerdoublefit{2});
                    
                    B1 = str2num(answerdoublefit{7});
                    
                    %Fit with a sum of two Lorentzians.
                    [f,~] = fit(xdata,ydata,Lorentz2,'StartPoint',[A2, A3, gamma2, gamma3, a2, a3, B1],'Weights',max(yerr)./yerr);
                    
                    % Plot the fit again so the user can look at it.
                    clf %Clear the current figure window
                    hold on
                    h0 = errorbar(xdata,ydata,yerr,'b.','Capsize',0.1); %Plot in the same figure window, overwriting the old one
                    h1 = plot(f);
                    
                    % Plot the guess as well to aid the user in fitting
                    h2 = plot(xdata, Lorentz2(A2,A3,gamma2,gamma3,a2,a3,B1,xdata), 'g');
                    legend([h0,h1,h2],'data', 'fit', 'guess')
                    
                    % Ask again if the fit is good.
                    answerfit = questdlg('Is this fit okay?','Manual fitting');
                    
                    switch answerfit
                        case 'No'
                            manualbadfit = 1;
                        case 'Cancel'
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case ''
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case 'Yes'
                            scantype(i) = 2; %Keep track of number of peaks for analysis later on
                            
                            % Get rid of the guess on the plot
                            clf
                            hold on
                            h0 = errorbar(xdata,ydata,yerr,'b.','Capsize',0.1);
                            h1 = plot(f);
                            legend([h0,h1],'data','fit')
                            
                            manualbadfit = 0; %Exit the while loop
                    end
                end
                j = j+1; %Index the while loop
            end
        end
        
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
            fitvalues(i).areas = [ values(1) intervals(1,1) intervals(2,1) ];
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
    
    % Provide feedback on how the spectra were fit.
    if feedbackglobal == 0
        fprintf(1, '\nAll data fit with a single Lorentzian.\n');
    elseif feedbackglobal == 1
        fprintf(1, '\nAC data fit with a sum of two Lorentzians.\n');
    end
    
     %% Plot all RPLE Spectra
    figure('Name', 'All Scans', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    cmap = jet(N);
    
    % Plot raw data.
    for i = 1:N
        plot(data(i).x,data(i).y,'Color',cmap(i,:))
    end
    
    legend(data(1:end).filename) %Legend labels spectra with filenames
    
    xlabel('Frequency (GHz)');
    ylabel('Intensity (arb. units)');
    hold off
    
    %% Make structure of arrays for fit parameters.
    % The fitvalues structure can tell you all things about one scan, now
    % we create a plotarrays structure that will tell you one thing about
    % all scans.
    
    plotarrays = struct('w0s',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1)),...
            'linewidths',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1)),...
            'heights',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1)),...
            'areas',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1)),...
            'B',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1)));
    
    % Preallocate space for arrays (increases code speed).
    indexref = zeros(1, N);
    indexAC = zeros(1, 2*N);
    
    % Reference
    plotarrays.w0s.ref = zeros(1, N);
    plotarrays.linewidths.ref = zeros(1, N);
    plotarrays.heights.ref = zeros(1, N);
    plotarrays.areas.ref = zeros(1, N);
    plotarrays.B.ref = zeros(1, N);
    
    % Reference errobars
    plotarrays.w0s.referr = zeros(1, N);
    plotarrays.linewidths.referr = zeros(1, N);
    plotarrays.heights.referr = zeros(1, N);
    plotarrays.areas.referr = zeros(1, N);
    plotarrays.B.referr = zeros(1, N);
    
    % AC
    plotarrays.w0s.AC = zeros(1, 2*N);
    plotarrays.linewidths.AC = zeros(1, 2*N);
    plotarrays.heights.AC = zeros(1, 2*N);
    plotarrays.areas.AC = zeros(1, 2*N);
    plotarrays.B.AC = zeros(1, 2*N);
    
    % AC error bars
    plotarrays.w0s.ACerr = zeros(1, 2*N);
    plotarrays.linewidths.ACerr = zeros(1, 2*N);
    plotarrays.heights.ACerr = zeros(1, 2*N);
    plotarrays.areas.ACerr = zeros(1, 2*N);
    plotarrays.B.ACerr = zeros(1, 2*N);
    
    % For loop to run through all scans
    for i = 1:N
        if scantype(i) == 0 %Reference scans
            indexref(i) = i;
            
            plotarrays.w0s.ref(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.ref(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.ref(i) = fitvalues(i).heights(1);
            plotarrays.areas.ref(i) = fitvalues(i).areas(1);
            plotarrays.B.ref(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.referr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2); %confint are symmetric so upper-lower to get error bar
            plotarrays.linewidths.referr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.referr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.referr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.referr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
        elseif scantype(i) == 1 %AC scan fit with 1 Lorentzian
            indexAC(i) = i;
            
            plotarrays.w0s.AC(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.AC(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.AC(i) = fitvalues(i).heights(1);
            plotarrays.areas.AC(i) = fitvalues(i).areas(1);
            plotarrays.B.AC(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.ACerr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2);
            plotarrays.linewidths.ACerr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.ACerr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.ACerr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.ACerr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
        elseif scantype(i) == 2 %AC scan fit with 2 Lorentzians
            % First peak of AC scans
            indexAC(i) = i;
            
            plotarrays.w0s.AC(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.AC(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.AC(i) = fitvalues(i).heights(1);
            plotarrays.areas.AC(i) = fitvalues(i).areas(1);
            plotarrays.B.AC(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.ACerr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2);
            plotarrays.linewidths.ACerr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.ACerr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.ACerr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.ACerr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
            % Second peak of AC scans
            indexAC(i + (N-1)) = i;
            
            plotarrays.w0s.AC(i + (N-1)) = fitvalues(i).w0s(4);
            plotarrays.linewidths.AC(i + (N-1)) = fitvalues(i).linewidths(4);
            plotarrays.heights.AC(i + (N-1)) = fitvalues(i).heights(4);
            plotarrays.areas.AC(i + (N-1)) = fitvalues(i).areas(4);
            plotarrays.B.AC(i + (N-1)) = plotarrays.B.AC(i);
            
            plotarrays.w0s.ACerr(i + (N-1)) = fitvalues(i).w0s(6)-fitvalues(i).w0s(5);
            plotarrays.linewidths.ACerr(i + (N-1)) = fitvalues(i).linewidths(6)-fitvalues(i).linewidths(5);
            plotarrays.heights.ACerr(i + (N-1)) = fitvalues(i).heights(6)-fitvalues(i).heights(5);
            plotarrays.areas.ACerr(i + (N-1)) = fitvalues(i).areas(6)-fitvalues(i).areas(5);
            plotarrays.B.ACerr(i + (N-1)) = plotarrays.B.ACerr(i);
        end
    end
    
    % Get rid of values that were unassigned.
    for i = N:-1:1 %Reference arrays
        if indexref(i) == 0
            indexref(i)=[];
        end
        if plotarrays.w0s.ref(i) == 0
            plotarrays.w0s.ref(i)=[];
        end
        if plotarrays.linewidths.ref(i) == 0
            plotarrays.linewidths.ref(i)=[];
        end
        if plotarrays.heights.ref(i) == 0
            plotarrays.heights.ref(i)=[]; 
        end
        if plotarrays.areas.ref(i) == 0
            plotarrays.areas.ref(i)=[]; 
        end
        if plotarrays.B.ref(i) == 0
            plotarrays.B.ref(i)=[]; 
        end
        if plotarrays.w0s.referr(i) == 0
            plotarrays.w0s.referr(i)=[];
        end
        if plotarrays.linewidths.referr(i) == 0
            plotarrays.linewidths.referr(i)=[];
        end
        if plotarrays.heights.referr(i) == 0
            plotarrays.heights.referr(i)=[];
        end
        if plotarrays.areas.referr(i) == 0
            plotarrays.areas.referr(i)=[];
        end
        if plotarrays.B.referr(i) == 0
            plotarrays.B.referr(i)=[];
        end
    end
    
    for i = 2*N:-1:1 %AC arrays
        if indexAC(i) == 0
            indexAC(i)=[];
        end
        if plotarrays.w0s.AC(i) == 0
            plotarrays.w0s.AC(i)=[];
        end
        if plotarrays.linewidths.AC(i) == 0
            plotarrays.linewidths.AC(i)=[];
        end
        if plotarrays.heights.AC(i) == 0
            plotarrays.heights.AC(i)=[];
        end
        if plotarrays.areas.AC(i) == 0
            plotarrays.areas.AC(i)=[];
        end
        if plotarrays.w0s.ACerr(i) == 0
            plotarrays.w0s.ACerr(i)=[];
        end
        if plotarrays.linewidths.ACerr(i) == 0
            plotarrays.linewidths.ACerr(i)=[];
        end
        if plotarrays.heights.ACerr(i) == 0
            plotarrays.heights.ACerr(i)=[];
        end
        if plotarrays.areas.ACerr(i) == 0
            plotarrays.areas.ACerr(i)=[];
        end
        if plotarrays.B.AC(i) == 0
            plotarrays.B.AC(i)=[];
        end
        if plotarrays.B.ACerr(i) == 0
            plotarrays.B.ACerr(i)=[];
        end
    end
    
    % Check that the reference arrays are correct length.
    if length(plotarrays.linewidths.ref) ~= Nref
        cprintf('err', '\nERROR: The length of the plotarrays.linewidths.ref array does not match the length of Nref.\n');
        beep; return
    end

    % Check that the AC arrays are correct length.
    if length(plotarrays.linewidths.AC) ~= length(indexAC)
        cprintf('err', '\nERROR: The length of the plotarrays.linewidths.AC array does not match the length of indexAC.\n');
        beep; return
    end
    
    %% Plot w0s vs. Scan Index.
    
    % Deal with the reference scans.
    if Nref == NAC
        smallpwr = 1;
        % Reference to each scan directly before (if magnitude of shift is small).
        refmean = zeros(1, 2*N);
        
        for i = 1:N %iterate through all AC scans
            if scantype(i) == 1
                refmean(i) = fitvalues(i-1).w0s(1);
            elseif scantype(i) == 2
                % First peak of AC scans
                refmean(i) = fitvalues(i-1).w0s(1);
                
                % Second peak of AC scans
                refmean(i + (N-1)) = fitvalues(i-1).w0s(1);
            end
        end
        
        % Clean array of values left unassigned.
        for i = 2*N:-1:1
            if refmean(i) == 0
                refmean(i)=[];
            end
        end
    else
        smallpwr = 0;
        % Proceed as usual with the average of all reference scans.
        refmean = mean(plotarrays.w0s.ref);
    end
    
    % Print the refmean to command line.
    fprintf(1, ['\nThe mean of reference scans is ' num2str(refmean) '.']);
    
    % Make a plot of w0s, with error bars for reference.
    figure('Name', 'w0s vs. Scan Index', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    errorbar(indexref,(plotarrays.w0s.ref-mean(refmean)),plotarrays.w0s.referr,'b.','markersize',15,'Capsize',0.1)
    errorbar(indexAC,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1)
    hline(0,'k:')
    
    legend('Reference', 'AC');
    
    xlabel('Scan Index');
    ylabel('w0s (GHz)');
    hold off
    
    %% User input.
    % Ask the user if they would like to input an array of power values, AC
    % Stark elevation angle, etc. to plot the w0s against.
    answer1 = questdlg('Would you like to input data to plot w0s against?', 'User data entry');
    
    switch answer1
        case 'No'
            userdataentered = 0;
        case 'Cancel'
            cprintf('err', '\nERROR: Please select either yes or no when asked to enter data.\n');
            beep; return
        case '' %User closed the dialog box
            cprintf('err', '\nERROR: Please select either yes or no when asked to enter data.\n');
            beep; return
        case 'Yes' %The user would like to enter an array of values
            userdataentered = 1; %Set for later
            
            % Settings for inputdlg() to ask user to enter data.
            prompt = {[ 'Enter array of length ' num2str(NAC) ':' ],'Data label:'};
            title = 'Data entry';
            dims = [1 100+NAC; 1 30];
            opts.Resize = 'on';
            opts.WindowStyle = 'normal';
            
            % If the program has been ran before use the previously entered
            % user data as the default input.
            if userenteredbefore == 0
                definput = {'[value1 value2 ... valueN]','AC Stark Laser Power (mW)'};
            elseif userenteredbefore == 1 && length(userdata) >= NAC
                definput = {num2str(userdata(1:NAC)),userdatalabel};
            elseif userenteredbefore == 1 && length(userdata) < NAC
                definput = {num2str(userdata),userdatalabel};
            end
            
            % Ask the user to enter the desired data.
            answer2 = inputdlg(prompt,title,dims,definput,opts);
            
            if isempty(answer2) %User pressed cancel
                cprintf('err', '\nERROR: Please enter the data when prompted.\n');
                beep; return
            else
                userdata = str2num(answer2{1});
                userdatalabel = answer2{2};
            end
    end
    
    if userdataentered == 1
        % Check that the user entered the approptiate number of values.
        while (length(userdata) ~= NAC) && (userdataentered == 1) %Keep prompting the user to enter data until they enter right amount
            % Tell the user they have entered the incorrect amount of data.
            opts.Default = 'OK';
            opts.Interpreter = 'none';
            answerA = questdlg('You have entered the incorrect number of values, please try again', 'ERROR',...
                'OK', 'No I give up', opts);

            switch answerA
                case 'No I give up' %The user has given up on entering data
                    userdataentered = 0; %This will cause the program to exit the while loop
                case '' %User closed the dialog box
                    cprintf('err', '\nERROR: That was rude. Please select one of the given options.\n');
                    beep; return
                case 'OK'
                    % Give the user another chance to enter the data.
                    prompt = {[ 'Enter array of length ' num2str(NAC) ': (The array previously entered was length '...
                     num2str(length(userdata)) ')' ],'Data label:'};
                    title = 'Data entry try again';
                    dims = [1 100+NAC; 1 30];
                    definput = {answer2{1},answer2{2}};
                    opts.Resize = 'on';
                    opts.WindowStyle = 'normal';
                    answer2 = inputdlg(prompt,title,dims,definput,opts);

                    if isempty(answer2) %User pressed cancel
                        cprintf('err', '\nERROR: Please enter the data when prompted.\n');
                        beep; return
                    else
                        userdata = str2num(answer2{1});
                        userdatalabel = answer2{2};
                    end
            end
        end
        % Perform some final cleanup on the userdata.
        
        % Double the length of userdata in case of double AC peaks.
        userdata = [ userdata userdata ];
        
        % Modify userdata based on if scans were fit with single or double
        % lorentzian.
        for i = N:-1:1
            if scantype(i) == 1
                ACscannum = sum(nAC(1:i));
                userdata(NAC + ACscannum) = [];
            end
        end
        
        % Check that userdata is now the correct length.
        if length(userdata) ~= length(indexAC)
            cprintf('err', '\nERROR: The length of the userdata array does not match the length of indexAC.\n');
            beep; return
        end
        
        % Equally space the reference scans throughout the userdata plot.
        step = range(userdata)/length(plotarrays.w0s.ref);
        userdataref = (min(userdata)+step/2):step:(max(userdata)-step/2);
    end
    
    %% Plot (w0s,linewidths,heights,areas) vs. User Data.
    
    if userdataentered == 1
        
        % w0s vs. userdata
        figure('Name', ['w0s vs. ' userdatalabel], 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        if smallpwr == 0
            h0 = errorbar(userdataref,(plotarrays.w0s.ref-refmean),plotarrays.w0s.referr,'b.','markersize',15,'Capsize',0.1);
            h1 = errorbar(userdata,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1);
        elseif smallpwr == 1
            h1 = errorbar(userdata,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1);
            h0 = hline(0,'k:');
        end
        hline(0,'k:')
        
        legend([h0 h1],'Reference', 'AC');
        
        xlabel(userdatalabel);
        ylabel('w0s (GHz)');
        hold off
        
        
        % Make another plot that includes the FWHM(gamma) values of each
        % peak plotted as a shaded area.
        figure('Name', ['w0s vs. ' userdatalabel ' with FWHM'], 'WindowStyle', 'docked', 'numbertitle', 'off');
        
        % Create arrays for +\- gamma/2. These are then plotted as lines,
        % and the are between them is shaded.
        if smallpwr == 0
            upperref = (plotarrays.w0s.ref-refmean)+plotarrays.linewidths.ref/2;
            lowerref = (plotarrays.w0s.ref-refmean)-plotarrays.linewidths.ref/2;
            upperAC = (plotarrays.w0s.AC-refmean)+plotarrays.linewidths.AC/2;
            lowerAC = (plotarrays.w0s.AC-refmean)-plotarrays.linewidths.AC/2;
        elseif smallpwr == 1
            upperAC = (plotarrays.w0s.AC-refmean)+plotarrays.linewidths.AC/2;
            lowerAC = (plotarrays.w0s.AC-refmean)-plotarrays.linewidths.AC/2;
        end
        
        % Calculate the number of AC scans fit with a single Lorentzian.
        nsinglepks = 2*NAC - length(plotarrays.w0s.AC);
        
        % Create the shaded regions for linewidths.
        hold on
        if smallpwr == 0
            h0 = fill([userdataref fliplr(userdataref)], [upperref fliplr(lowerref)], [0.93 0.93 0.93], 'linestyle', 'none');
            fill([userdata(1:NAC) fliplr(userdata(1:NAC))],...
                [upperAC(1:NAC) fliplr(lowerAC(1:NAC))], [0.93 0.93 0.93], 'linestyle', 'none');
            fill([userdata((end/2)+(nsinglepks/2)+1:end) fliplr(userdata((end/2)+(nsinglepks/2)+1:end))],...
                [upperAC((end/2)+(nsinglepks/2)+1:end) fliplr(lowerAC((end/2)+(nsinglepks/2)+1:end))], [0.93 0.93 0.93], 'linestyle', 'none')
        elseif smallpwr == 1
            h0 = fill([userdata(1:NAC) fliplr(userdata(1:NAC))],...
                [upperAC(1:NAC) fliplr(lowerAC(1:NAC))], [0.93 0.93 0.93], 'linestyle', 'none');
            fill([userdata((end/2)+(nsinglepks/2)+1:end) fliplr(userdata((end/2)+(nsinglepks/2)+1:end))],...
                [upperAC((end/2)+(nsinglepks/2)+1:end) fliplr(lowerAC((end/2)+(nsinglepks/2)+1:end))], [0.93 0.93 0.93], 'linestyle', 'none')
        end
        
        % Add the w0s to the plot.
        hold all
        if smallpwr == 0
            h1 = errorbar(userdataref,(plotarrays.w0s.ref-refmean),plotarrays.w0s.referr,'b.','markersize',15,'Capsize',0.1);
            h2 = errorbar(userdata,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1);
            hline(0,'k:');
        elseif smallpwr == 1
            h1 = hline(0,'k:');
            h2 = errorbar(userdata,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1);
        end
        legend([h0 h1 h2], '\gamma  (FWHM)', 'Reference', 'AC');
        
        xlabel(userdatalabel);
        ylabel('w0s (GHz)');
        hold off
        
        
        % linewidths vs. userdata
        figure('Name', ['Linewidths vs. ' userdatalabel], 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        errorbar(userdataref,plotarrays.linewidths.ref,plotarrays.linewidths.referr,'b.','markersize',15,'Capsize',0.1)
        errorbar(userdata,plotarrays.linewidths.AC,plotarrays.linewidths.ACerr,'r.','markersize',15,'Capsize',0.1)
        
        legend('Reference', 'AC');
        
        xlabel(userdatalabel);
        ylabel('Linewidth (GHz)');
        hold off
        
        
        % heights-background vs. userdata
        figure('Name', ['Peak Heights - Background vs. ' userdatalabel], 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        errorbar(userdataref,plotarrays.heights.ref-plotarrays.B.ref,plotarrays.heights.referr,'b.','markersize',15,'Capsize',0.1)
        errorbar(userdata,plotarrays.heights.AC-plotarrays.B.AC,plotarrays.heights.ACerr,'r.','markersize',15,'Capsize',0.1)
        
        legend('Reference', 'AC');
        
        xlabel(userdatalabel);
        ylabel('Peak Height (Int)');
        hold off
        
        
        % areas vs. userdata
        figure('Name', ['Peak Areas vs. ' userdatalabel], 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        errorbar(userdataref,plotarrays.areas.ref,plotarrays.areas.referr,'b.','markersize',15,'Capsize',0.1)
        errorbar(userdata,plotarrays.areas.AC,plotarrays.areas.ACerr,'r.','markersize',15,'Capsize',0.1)
        
        legend('Reference', 'AC');
        
        xlabel(userdatalabel);
        ylabel('Peak Area (Int*GHz)');
        hold off
    end
    
    %% Report average/max shift to the command line.
    
    % Average shift if the user has not entered data.
    if feedbackglobal == 0 && userdataentered == 0
        ACmean = mean(plotarrays.w0s.AC);
        fprintf(1, ['\nThe average shift is ', num2str(ACmean-refmean), ' GHz.\n']);
    elseif feedbackglobal == 1 && userdataentered == 0
        ACmean1 = mean(plotarrays.w0s.AC(1:length(plotarrays.w0s.AC)/2));
        ACmean2 = mean(plotarrays.w0s.AC(1+length(plotarrays.w0s.AC)/2:length(plotarrays.w0s.AC)));
        fprintf(1, ['\nThe average shift of peak 1 is ', num2str(ACmean1-refmean), ' GHz.']);
        fprintf(1, ['\nThe average shift of peak 2 is ', num2str(ACmean2-refmean), ' GHz.\n']);
    end
    
    % Max shift if user has entered data.
    if userdataentered == 1
       maxshift =  sign(max(plotarrays.w0s.AC-refmean))*max(abs(plotarrays.w0s.AC-refmean));
       fprintf(1, ['\nThe maximum shift is ', num2str(maxshift), ' GHz.\n']);
    end
    
    %% Save the figure window and the workspace.
    
    % Array with the figure handles.
    Nplots = N + 5*userdataentered + 2;
    h = 1:Nplots;
    
    %Save the figure.
    savefig(h, [ path '\' folder{1}, ' Plots', '.fig' ]);
    
    %Save specific variables from the workspace.
    if feedbackglobal == 0 && userdataentered == 0
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'ACmean', 'data', 'fitvalues', 'indexAC', 'indexref', 'refmean',...
        'NAC', 'Nref', 'scantype', 'usertolerance', 'plotarrays');
    elseif feedbackglobal == 0 && userdataentered == 1
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'indexAC', 'indexref', 'refmean',...
        'userdata', 'userdatalabel', 'userdataref', 'NAC', 'Nref', 'scantype', 'maxshift', 'usertolerance','plotarrays');
    elseif feedbackglobal == 1 && userdataentered == 0
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'ACmean1', 'ACmean2', 'data', 'fitvalues', 'indexAC', 'indexref', 'refmean',...
        'NAC', 'Nref', 'scantype', 'usertolerance', 'plotarrays');
    elseif feedbackglobal == 1 && userdataentered == 1
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'indexAC', 'indexref', 'refmean',...
        'userdata', 'userdatalabel', 'userdataref', 'NAC', 'Nref', 'scantype', 'maxshift', 'usertolerance', 'plotarrays');
    end
    
    fprintf(1, '\nData and figures saved!\n');
    
end