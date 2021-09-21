% ******************************
Velocity=xlsread('FUDS');
V=Velocity(:,2)*1.6;
N=length(V); % No. of Samples
mass = [180 600 1500]; % 2W, 3W & 4W 
area = [0.6 1.6 2.3]; % Frontal area in square metres
Cd = [1 0.45 0.32] ; % Drag coefficient
Gratio = [2/0.28 5/0.2 8/0.3]; % Gearing ratio, = G/r
Pac=[50 200 500]; % Average power of accessories in Watt

G_eff = 0.95; % Transmission efficiency
eff_mot=0.80;%Motor Efficiency
Regen_ratio = 1; %This sets the proportion of the

%Initialization of Variables
D_end = zeros(1,100);
Pbat=zeros(1,N);
Pregen=zeros(1,N);
D=zeros(1,N); % Record of distance travelled in km.
Ebat=zeros(1,N);
PTE=zeros(1,N);
CY=1;

i_type=3;% 1 for 2W, 2 for 3W and 3 for 4W

for C=2:N
    accel=(V(C) - V(C-1))/3.6;%Acceleration
    D(C) = D(C-1) + (V(C)/3.6/1000);%Distance
    
    Fad = 0.5 * 1.25 * area(i_type) * Cd(i_type) * (V(C)/3.6)^2;% Aero Drag
    Frr=0.015 * mass(i_type) * 9.8; % Rolling Resistance
    Fhc = 0; % Hill Climb Resistance
    Fla = mass(i_type) * accel; % Linear Acceleration
    Pte = (Frr + Fad + Fhc + Fla)*(V(C)/3.6);%Traction Power
    PTE(C)=Pte; % Traction Power
    omega = Gratio(i_type) * (V(C)/3.6);%Motor rotational speed

    
    if omega == 0 % Stationary
        Pte=0;
        Pmot_in=0; % No power into motor
        Torque=0; % No Torque
        
    elseif omega > 0 % Moving        
        if Pte>=0 % Crusing or Accelerating
            Pmot_out=Pte/G_eff; % Motor power > shaft power
            Pmot_in = Pmot_out/eff_mot; % Motoring Mode
        elseif Pte<0 % Braking 
            Pmot_out= Regen_ratio * Pte * G_eff; % Motor power < shaft power
            Pmot_in = Pmot_out * eff_mot; % Generating Mode
        end        
    Torque=Pmot_out/omega; % Torque: P=T*omega
    end
    
    Pbat(C) = Pmot_in + Pac(i_type);   % Power required from Battery
    Ebat(C)=Ebat(C-1)+Pbat(C); % Energy Required from Battery
    
    
    TIME(C)=C; % Time
    VEL(C)=V(C); %Velocity (Drive Cycle)
    ACC(C)=accel; %Acceleration (Drive Cycle)
    P_MOT(C)=Pmot_in/1000; % Power from Motor
    OMEGA(C)=omega; % Motor Speed in rpm=Wheel Speed *G
    T_MOT(C)=Torque; % Motor Torque
    DIST(C)=D(C); % Current from Battery
    ENERGY_B(C)=Ebat(C)/(3600*1000);
    
    Omegamax=max(OMEGA);
    Vmax=max(VEL);
    Accmax=max(ACC);
    Tmax=max(T_MOT);
    Pmax=max(P_MOT);
    Pavg=mean(P_MOT);    
end
figure(1)
subplot(3,1,1);
plot(TIME,VEL);
xlabel('Time(sec)');
ylabel('Vehicle Speed (km/hr)');
str=['Max Velocity :',num2str(Vmax),' km/hr'];
title(str);
subplot(3,1,2);
plot(TIME,ACC);
xlabel('Time(sec)');
ylabel('Vehicle Acceleration(m/s2)');
str=['Max Acc :',num2str(Accmax),' m/s2'];
title(str);
subplot(3,1,3);
plot(TIME,DIST);
xlabel('Time(sec)');
ylabel('Distance covered (km)');
str=['Distance travelled :',num2str(DIST(:,end)),' km'];
title(str);
figure(2)
subplot(3,1,1);
plot(TIME,T_MOT);
xlabel('Time(sec)');
ylabel('Motor Torque Req (Nm)');
str=['Motor(Peak Torque): ',num2str(Tmax),' Nm'];
title(str);
subplot(3,1,2);
plot(TIME,OMEGA);
xlabel('Time(sec)');
ylabel('Motor Speed Req (rad/sec)');
str=['Motor(Peak Speed): ',num2str(Omegamax),' rad/sec'];
title(str);
subplot(3,1,3);
plot(TIME,P_MOT);
xlabel('Time(sec)');
ylabel('Motor Power Req (kW)');
str=['Motor(Peak Power): ',num2str(Pmax),' kW & (Avg Power): ',num2str(Pavg),' kW'];
title(str);
figure(3)
plot(TIME,ENERGY_B);
xlabel('Time(sec)');
ylabel('Battery Energy Req (kWh)');
str=['Battery Energy Req : ',num2str(ENERGY_B(:,end)),' kWh'];
title(str);