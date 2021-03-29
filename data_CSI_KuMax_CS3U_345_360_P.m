clear all
clc

% http://www.ecomp.uefs.br/public_html/formularios.php
% Reunião 14/03/2021

%% Information from the CSI_KuMax_CS3U-345-360-P solar array datasheet
Isc = 9.67; % Corrente de Curto Circuito na STC [A]
Voc = 47.00; % Tensão de Circuito Aberto na STC [V] 
Imp = 9.10; %Array current @maximum power point [A]
Vmp = 39.6; %Array maximum output peak power [W]
Kv = -0.29*10^-2; %Voltage/ temperature coefficient [V/K]
Ki = 0.05*10^-2; %Current/ temperature coefficient [A/K]
Ns = 144; %Number of series cells
G =  1000; %Nominal radiance [W/m^2] @ 25ºC
Tn = 25 + 273.15; %Nominal operating temperature [K]
n = 1.3; #The ideality factor of diode
SilAmorfoV = 0.8426;
SilAmorfoI = 0.9411;

T = Tn;

%% Constants
K = 1.3806503e-23; %Boltzmann [J/K]
q = 1.60217646e-19; %Electron charge [C]
%% a = 1.0; %Ideal diod

% Calcular o Ipv
% Calcular o Io (Com base no Irs)
% Estimar valor para Rp e Rs
% Gerar um vetor de V
%Calcular o I através de Newton-Raphson, I = f(I,V) Output current , passa
%tudo pro outro lado e iguala a zero

%% Iph

Iph = (Isc + Ki * (T - 298)) * (G/1000);

%% Io
Ego = 1.1; #band gap energy of the semiconductor
Irs = Isc/(e^((q * Voc)/(n * Ns * K * T)) - 1);
Io = Irs * (T/Tn)^3 * exp((q * Ego * (1/Tn - 1/T))/(n * K));

%% Vetor de V

%Voltage points in each iteraction
nv = 100; %Defines how many points are used for obtaining the IxV curve

V = 0:Voc/nv:Voc;
%I = zeros(1,size(V,2));    % Current vector
Rs = 0.221; %Serie resistence Ohm
Rsh = 405.415; %Shunt resistence Ohm
%% I

#I = Iph - Io * (exp((q*(V + I * Rs))/(n * K * Ns * T))-1) - ((V + I * Rs)/Rsh);
I = @(v)Iph - Io .* (exp((q.*(v + I .* Rs))./(n .* K .* Ns .* T))-1) - ((v + I .* Rs)./Rsh);
potencia = @(v) v .*I(v); # f(v) = v * I(v)
dIdf = derivada(potencia); #primeira derivada
dI2df = derivada(dIdf); #segunda derivada - Deu algum bug nesta segunda derivada.
#Possivelmente, as casas decimais relevantes ultrapassaram as casas decimais do epsilon da m�quina
#newton(dIdf, dI2df, 30 , 0.001, 100);
secante(dIdf,30,47,0.0000001,100);

%% SAÍDAS GRÁFICAS

figure(1)
plot(V, dIdf(V))
#plot(V, I(V))
#plot(V,I)
xlabel('Tensão Terminal, [V]')
ylabel('Potência, [W] (Derivada)')

figure(2)
#plot(V,dI2df(V))
plot(V,potencia(V))
#plot(V,I(V))
xlabel('Tensão Terminal, [V]')
ylabel('Potência, [W]')

figure(3)
#plot(V,dI2df(V))
#plot(V,potencia(V))
plot(V,I(V))
xlabel('Tensão Terminal, [V]')
ylabel('Corrente, [A]')


%Vtn = k * Tn / q;           %Thermal junction voltage (nominal) 
%Vt  = k * T  / q;           %Thermal junction voltage (current temperature)


