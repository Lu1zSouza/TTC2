%% PV modeling - Method 1 - COMPLETE MODEL WITH Rs and Rp

% Author: Marcelo Gradella Villalva

% University of Campinas - UNICAMP - Brazil - 2015
% mvillalva@gmail.com

% If you like my work please cite my papers:

% M. G. Villalva, J. R. Gazoli, E. Ruppert F., 
% "Comprehensive approach to modeling and simulation 
% of photovoltaic arrays", IEEE Transactions on 
% Power Electronics, 2009, % vol. 25, no. 5, 
% pp. 1198--1208, ISSN 0885-8993
 
% M. G. Villalva, J. R. Gazoli, E. Ruppert F.
% "Modeling and circuit-based simulation of photovoltaica arrays"
% Brazilian Journal of Power Electronics, 2009
% vol. 14, no. 1, pp. 35--45, ISSN 1414-8862

% Visit for updates: http://sites.google.com/site/mvillalva/pvmodel


%% Load PV datasheet information

clear all

 %data_KC85T
 %data_KD135SX_UPU
 %data_MSX60
 %data_amorphous
 %data_MSX60_single_cell
 %data_Q6LM
 data_M2453BB
 

%% PV model calculation 

% Warning: do not change any values here unless you know what you are
% doing. You will have the opportunity to use any value of temperature or
% irradiance when you evaluate the model using the PV_model_eval program.

Gn = 1000;               %Nominal irradiance [W/m^2] @ 25oC
Tn = 25 + 273.15;        %Nominal operating temperature [K]

% Parameters used to control how the algorithm will work

%Increment of Rs
Rsinc = 0.001;

% Rsinc controls the precision and the algorithm speed. Small values are
% better for improved precision and helps the program to converge

%Maximum tolerable power error
tol = 0.0001; %Defines the model precision

%Number of "a" values to try
%namax = 100; %This will be used in a future version of the program

%Increment of "a"
%ainc = 0.01; %This will be used in a future version of the program

%Voltage points in each iteraction
nv = 100; %Defines how many points are used for obtaining the IxV curve

%Maximum number of iteractions for each value of "a"
%Avoids  program stall in case of non convergence
%This will be useful in a future version of the program
nimax = 500000; 

%used for debugging
plott = 0; %1 = Enables plotting during algorithm execution
           %0 = Disables plotting (better)


%% PROGRAM STARTS HERE

% Modeling algorithm - here we are obtaining the PV model parameters

% Reference values of Rs and Rp  
% These values are not used in the modeling proces but they will 
% be displayed at the end. 
Rs_max = (Vocn - Vmp)/ Imp;
Rp_min = Vmp/(Iscn-Imp) - Rs_max;

% Initial guesses of Rp and Rs
Rs = 0;
Rp = Rp_min;

% The model is adjusted at the nominal condition
T = Tn;
G = Gn; 

k = 1.3806503e-23;   %Boltzmann [J/K]
q = 1.60217646e-19;  %Electron charge [C]

Vtn = k * Tn / q;           %Thermal junction voltage (nominal) 
Vt  = k * T  / q;           %Thermal junction voltage (current temperature)

perror = Inf; %dummy value
ni = 0; %counter
a = 1; %Initial value of a

% Iterative loop executed untin Pmax,model = Pmax,experimental
while (perror>tol) && (Rp > 0) && (ni < nimax)

ni = ni + 1; 
    
%  Temperature and irradiation effect on the current
dT = T-Tn;
Ipvn = (Rs+Rp)/Rp * Iscn;         % Nominal light-generated current
Ipv = (Ipvn + Ki*dT) *G/Gn;       % Actual light-generated current 
Isc = (Iscn + Ki*dT) *G/Gn;       % Actual short-circuit current

Ion = (Ipv - Vocn/Rp)/(exp(Vocn/Vt/a/Ns)-1); 
Io = Ion;

% Increments Rs 
Rs  = Rs + Rsinc;  
Rp_ = Rp;

% Egap = 2.72370016e-19;  % Bandgap do silício amorfo em J (=1.7 eV)
  Egap = 1.8e-19;         % Bandgap do silício cristalino em J (=1.124 eV)

a = (Kv - Vocn/Tn) / ( Ns * Vtn * ( Ki/Ipvn - 3/Tn - Egap/(k*Tn^2) ) );
 
% Comments:

% In previous versions of this program I used to determine the value of
% "a" by trial and error, i.e. I chose some value of "a" before running
% the program. In order to obtain the parameters of the two-diode model
% you need to guess either the value of "a" or the value of Egap for the
% device you are modeling. In this version of the program I'm using the
% value of Egap of Si found in the litterature and calculating "a"
% during the loop execution. This worked nicely with the devices 
% I've tested. If you find any problems just uncomment the line bellow
% with any value of "a" you wish.

% a = 1; % Read the comments above.
% You can try other values of the ideality factor if your results are not 
% satisfactory or if you find convergence problems.

Rp = Vmp*(Vmp+Imp*Rs)/(Vmp*Ipv-Vmp*Io*exp((Vmp+Imp*Rs)/Vt/Ns/a)+Vmp*Io-Pmax_e);

% Solving the I-V equation for several (V,I) pairs 
clear V
clear I

V = 0:Vocn/nv:Vocn;        % Voltage vector
I = zeros(1,size(V,2));    % Current vector

for j = 1 : size(V,2) %Calculates for all voltage values 
    
% Solves g = I - f(I,V) = 0 with Newton-Raphson
  
g(j) = Ipv-Io*(exp((V(j)+I(j)*Rs)/Vt/Ns/a)-1)-(V(j)+I(j)*Rs)/Rp-I(j);
  
while (abs(g(j)) > 0.001)
      
g(j) = Ipv-Io*(exp((V(j)+I(j)*Rs)/Vt/Ns/a)-1)-(V(j)+I(j)*Rs)/Rp-I(j);
glin(j) = -Io*Rs/Vt/Ns/a*exp((V(j)+I(j)*Rs)/Vt/Ns/a)-Rs/Rp-1;
I_(j) = I(j) - g(j)/glin(j);
I(j) = I_(j);   
  
end  

end % for j = 1 : size(V,2)

if (plott) 
    
  %Plots the I-V and P-V curves

  %Current x Voltage
  figure(1) 
  grid on
  hold on 
  title('I-V curve - Adjusting Rs and Rp');
  xlabel('V [V]');
  ylabel('I [A]');
  xlim([0 Vocn]);
  ylim([0 Iscn]);
 
  %Plots I x V curve
  plot(V,I,'LineWidth',2,'Color','k') 
 
  %Plots the "remarkable points" on the I x V curve
  plot([0 Vmp Vocn],[Iscn Imp 0],'o','LineWidth',2,'MarkerSize',5,'Color','k') 

  %Power x Voltage
  figure(2)  
  grid on
  hold on  
  title('P-V curve - Adjusting peak power');
  xlabel('V [V]');
  ylabel('P [W]');
  xlim([0 Vocn])
  ylim([0 Vmp*Imp]);
    
end % if(plott)

  % Calculates power using the I-V equation
  P = (Ipv-Io*(exp((V+I.*Rs)/Vt/Ns/a)-1)-(V+I.*Rs)/Rp).*V;
  
  Pmax_m = max(P);
  
  perror = (Pmax_m-Pmax_e);
  
if (plott) 
  
  %Plots P x V curve
  plot(V,P,'LineWidth',2,'Color','k')
  
  %Plots the "remarkable points" on the power curve
  plot([0 Vmp Vocn],[0 Vmp*Imp 0],'o','LineWidth',2,'MarkerSize',5,'Color','k')
  
end % if (plott)
     
end % while (error>tol) 

if (Rp<0) Rp = Rp_ 
end

% PROGRAM ENDS HERE

%% Outputs
 
 % I-V curve
 figure(3) 
 grid on
 hold on 
 title('Adjusted I-V curve');
 xlabel('V [V]');
 ylabel('I [A]');
 xlim([0 max(V)*1.1]);
 ylim([0 max(I)*1.1]);
 plot(V,I,'LineWidth',2,'Color','k') %
 plot([0 Vmp Vocn ],[Iscn Imp 0 ],'o','LineWidth',2,'MarkerSize',5,'Color','k') 
 
 % P-V curve
 figure(4) 
 grid on
 hold on 
 title('Adjusted P-V curve');
 xlabel('V [V]');
 ylabel('P [W]');
 xlim([0 Vocn*1.1]);
 ylim([0 Vmp*Imp*1.1]);  
 plot(V,P,'LineWidth',2,'Color','k') %
 plot([0 Vmp Vocn ],[0 Pmax_e 0 ],'o','LineWidth',2,'MarkerSize',5,'Color','k') 
   
disp(sprintf('Method 1 - complete model\n'));
disp(sprintf(' Rp_min = %f',Rp_min));
disp(sprintf('     Rp = %f',Rp));
disp(sprintf(' Rs_max = %f',Rs_max));
disp(sprintf('     Rs = %f',Rs));
disp(sprintf('      a = %f',a));
disp(sprintf('      T = %f',T-273.15));
disp(sprintf('      G = %f',G));
disp(sprintf(' Pmax,m = %f  (model)',Pmax_m));
disp(sprintf(' Pmax,e = %f  (experimental)',Pmax_e));
disp(sprintf('    tol = %f',tol));
disp(sprintf('P_error = %f',perror));
disp(sprintf('    Ipv = %f',Ipv));
disp(sprintf('    Isc = %f',Isc));
disp(sprintf('    Ion = %g',Ion));
disp(sprintf('\n\n'));

%%% This feature is not available yet.
%disp('Copy and paste this line of code in your simulation model:');
%disp(sprintf('pvdata[14]={%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f};',Iscn,Vocn,Imp,Vmp,Kv,Ki,Ns,Gn,Tn,Rp,Rs,a,Ipvn,Ion));

