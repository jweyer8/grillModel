%%Heat Transfer Grilling Project

%clean up
clear all;
close all;
clc;

%%Contstant vars

%tofu constants 
a = 20; %accuracy/number of nodes
thickness = ((1)*2.54)/100; %thickness of the tofu [m]
L = thickness; %tofu length [m]
burnTemp = 370; % this is the boling temp of water and the temperature that the tofu could burn [K]
idealTemp = 320; %ideal iner temperature of tofu [K]
tofuStartTemp = 290; % starting temperature of tofu [K]

%air constants
infTemp = (180 - 32) * 5/9 + 273.15; %grill air temp [K]

%grill cotnstants
stoveHeatFlux = 5000; %heat flud from grill [W/m^2]

%other constants
delt = 7.5; %time step [s]
dely = thickness/a; % space between each node [m]
time = 0; %start time
downSide = 1; %side of the tofu that is facing the grill
flipCount = 0; %number of times tofu has been flipped

%%1D Implicite Finite Differnce Matrix 

A = zeros(a,a); %implicit finite differnce matrix
T = zeros(a,1); %temperatures of the nodes for specific time
C = zeros(a,1);

%set temp of nodes equal to the starting temp of the tofu
for i = 1:a
   T(i,1) = tofuStartTemp;
end


%%save centerline, and outer temps through time

center = [];
outer1 = [];
outer2 = [];
timeArray = [];



%%set up for animation
hh = figure;
axis tight manual 
filename = 'C:\Users\jweye\Documents\Fluids\grillProject\tempAnimation.gif';



%%Marching through time 

%while time is less than given time keep looping through values
while(T(round(a/2)) < idealTemp)
    
    %recaulate the finite differnce matrix with new temperatures
    [A,F,B,k] = createFinite(T,infTemp,delt,dely,L);
    
    %set C matrix equal to the previous node temps
    for i = 1:a
        C(i) = T(i);
    end
    
    %create new C matrix based off previous temp and other constants
    if(downSide == a)
        C(a) = C(a) + 2*F*stoveHeatFlux*dely/k; %add constant heat flux
        C(1) = C(1) + 2*F*B*infTemp; %add convective term
    else
        C(1) = C(1) + 2*F*stoveHeatFlux*dely/k;
        C(a) = C(a) + 2*F*B*infTemp; 
    end
    
    %calculate new node temps
    T = inv(A)*C;
    
    %flip tofu
    if(T(downSide) > burnTemp && flipCount < 2) 
        A = flipud(A);
        A = fliplr(A);
        flipCount = flipCount + 1;
        fprintf('\nTofu fliped at %.3f minute mark\n',time/60);
        
        if(downSide == a)
            downSide = 1;
        else
            downSide = a;
        end
    end
    
    %save temps to respective arrays so we can plot later
    center = [center, {T(round(length(T)/2))}];
    outer1 = [outer1, {T(1)}];
    outer2 = [outer2, {T(a)}];
    timeArray = [timeArray, time];
   
    %animation of temps
    TT = [T,T,T,T,T];
    contourf(TT,a,'LineStyle','none')
    colormap hot
    axis off;
    drawnow
    
    %save frames
    frame = getframe(hh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if time == 0
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    %increase time based off the given time step
    time = time + delt;
end

%%Plotting

%output some data about the cooking of tofu
min = time/60;
fprintf('The number of times tofu was fliped is: %.1f\n',flipCount);
fprintf('\nThe time it took to cook the tofu is: %.1f minutes \n\n',min);
fprintf('The temperature of the nodes ate the end of cooking time are:\n\n');
T

%graph temperature vs time for centerline and outside nodes
figure;
hold on;
plot(timeArray,cell2mat(center),'linewidth',2);
plot(timeArray,cell2mat(outer1),'linewidth',2);
plot(timeArray,cell2mat(outer2),'linewidth',2);
ylim([tofuStartTemp-10 450]);
legend('centerline temperature','bottom (before flips)','top (before flips)','location','northwest');
title('Temperature vs Time');
ylabel('Temperature [K]');
xlabel('Time [s]');
hold off
 

%%finite differnce based on varying temperatures of tofu
 
function [A,F,B,k] = createFinite(T,infTemp,delt,dely,L)
    
    a = 20; %acuracy/number of nodes
    g = 9.81; %gravity [m/s^2]
    moi = .5; %moisture content
    range = 25;
    Tavg = 27; %average temp value of tofu [C]
    kuAir = 15.89*10^-6; %kenematic viscosity of air [m^2/s]
    diffuseAir = 22.5*10^-6; %thermal diffusivity of air [m^2/s]
              
    
    %tofu thermal properties based off average temp of tofu 
    k = (.2112 + (8.943*10^-4)*moi*Tavg + .3077*(moi^2)) %conduction coefienct of tofu based on moisture content and average temp [W/mK]
    diffuse = (.0816 - .05682*moi + .1164*(moi^2) + (6.866*10^-4)*(moi^2)*Tavg - (5.17*10^-6)*(moi^2)*(Tavg^2))*10^-6%thermal diffusivity of tofu based on moisture content and average temp

    
    %calculated constants
    filmTemp = (T(a) + infTemp)/2; %film temperature, this will change as tofu is cooked however it will have neglible change on the convection coeffient [K]
    beta = 1/(filmTemp); %Boussinesq Number with ideal gas approx
    RA = abs((g*beta*(T(a) - infTemp)*L^3)/(kuAir*diffuseAir)) %Raligh number

    if(RA < 10^7)
        N  = .54 * RA^(1/4); %nusselt below critical
    else
        N = .15 * RA^(1/3); %nusselt above critical
    end
    
    h = 10;% N*kAir/L %free convective flow coefient [W/m^2K]
    F = (diffuse*delt)/dely^2 %Fourier Number
    B = (dely*h)/k %Biot Number


    %create the 1D implicite finite difference array (time independent)
    for i = 1:a
        %center nodes (two conduction)
        if(i ~= a && i ~= 1)
            A(i,i) = 1 + 2*F;
            A(i,i+1) = -F;
            A(i,i-1) = -F;
        %node touching stove (one conduction & one constant heat flux)
        elseif(i == 1)
            A(i,i) = 1 + 2*F;
            A(i,i+1) = -2*F;
        %node touching air (one conduction & one convection)
        else
            A(i,i) = 1 + 2*F + 2*B*F;  
            A(i,i-1) = -2*F;
        end
    end
end





