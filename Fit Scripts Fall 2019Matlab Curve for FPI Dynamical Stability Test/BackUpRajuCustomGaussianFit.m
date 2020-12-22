%% Upload some data which will be used as the dependent and the independent 
%variable in the Matlab.
x = Data1(:,1); % x column of the data
y = Data1(:,2); % y column of the data

%% Rough initial guesses I dont know how to choose the rough initial guess
%a3 = ((max(x)- min(x))/10)^2;
%a2 = ((max(x) + min(x)))/2;
%a1 = max(y)*a3;
%a0 = [a1,a2,a3];  % Initial guess parameters

%% define the fitting function with fittype
% notice that it is quite human readable
% Matlab automatically treats x as independent variable
%f=fittype(@(A,x0,gamma, x) (A./(pi*gamma))*(1+((x-x0).^2)./gamma).^-1);
% f=fittype(@(A,x0,gamma,B,x)(A/pi)*(0.5*gamma)*(((x-x0).^2+0.25*(gamma^2)).^-1)+B);
 f=fittype(@(A,x0,gamma,B,x)(A/pi)*(0.5*gamma)*(((x-x0).^2+0.25*(gamma^2)).^-1)+B);


%% assign initial guessed parameters
% [A, x0, gamma] they are in the order of the appearance
% in the above fit function definition
%pin=[7.276 0.02309 0.0007 0.005]; % Initial guesses to fit the Gaussian function
pin=[7.276 0.02307 0.0020 0]; % Initial guesses to fit the Gaussian function

%% Finally, we are ready to fit our data
[fitobject] = fit (x,y, f, 'StartPoint', pin);

%% Let's see how well our fit follows the data
plot(fitobject, x,y,'--k', 'fit')
set(gca,'FontSize',24); % adjusting font size
xlabel('x');
ylabel('y');