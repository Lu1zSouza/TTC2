clear all
clc


%% Information from the KC85T solar array datasheet

Iscn = 5.34;             %Nominal short-circuit voltage [A] 
Vocn = 21.7;            %Nominal array open-circuit voltage [V] 
Imp = 5.02;              %Array current @ maximum power point [A] 
Vmp = 17.4;             %Array voltage @ maximum power point [V] 
Pmax_e = Vmp*Imp;       %Array maximum output peak power [W] 
Kv = -8.21e-2;            %Voltage/temperature coefficient [V/K] 
Ki = 2.12e-3;              %Current/temperature coefficient [A/K] 
Ns = 36;                %Nunber of series cells 
Gn = 1000;              %Nominal irradiance [W/m^2] @ 25oC
Tn = 25 + 273.15;       %Nominal operating temperature [K]


%% Constants

k = 1.3806503e-23;   %Boltzmann [J/K]
q = 1.60217646e-19;  %Electron charge [C]

a = 1.0;

%% Algorithm parameters

%Increment of Rs
Rsinc = 0.0001;

%Initial value of "a"
%a = 1.0; 

%Increment of "a"
%ainc = 0.01;

%Maximum tolerable power error
tol = 0.00001;

%Maximum number of iteractions for each value of "a"
nimax = 100000;

%Voltage points in each iteraction
nv = 2000;

%Number of "a" values to try
%namax = 50;




%% Experimental points collected from datasheet

% not available


