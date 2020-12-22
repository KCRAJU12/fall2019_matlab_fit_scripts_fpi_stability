%% The following script plots the second column of the extracted coefficients from the 
%% fitted parameter matrices which is named as C

[FileName,PathName] = uigetfile('*.mat', 'Open mat file','MultiSelect','on');
 data = cell(length(FileName),1) ;
 
 for i = 1:length(FileName) 
     file = load(fullfile(PathName,FileName{i}));
     data{i} = file ;
 end

 figure;
 hold on; % This waits for until the for loop stops. If we do not include this line in each iteration of the for loop
 % the previous value is overwritten. So this line is important to include.
 

 Maxvalue = cell(length(FileName),1); %  maximum value of scan range
 
 Minvalue = cell(length(FileName),1); %  minimum value of the scan range
 
 for k = 1: length(FileName)
     
     xdata = linspace(0,80,244);
     
     ydata = data{k,1}.C(1:244,2);
     
     plot(xdata,ydata); % The order of this line matters here 
     
     title('FPI dynamical stability test');
     
     xlabel('Time (s)');
     
     ylabel('Scan range (s)');
     
     Maxvalue{k} = max(data{k,1}.C(1:244,2));
    
     Minvalue{k} = min(data{k,1}.C(1:244,2));
 end
 
 
 hold off