%% This function tends to establish a relationship between coordinates in LVLH frame and ECI frame 
%% simulate PCO and GCO formation flying in ECI and LVLH frame
function ECI_LVLHcoordinatetransform
%satellite initial unit velocity vector and velocity magnitude
satrunspeed=300;
%NUMBER OF REVOLUTIONS 
rev_number=5;
mu = 398.6004418e12;  % Planetary gravitational constant for Earth, (mu = GMearth) (m^3/s^2)
%Number of Days since Jan 1, 2000
J2000_days=103752/24; % = 4321 on 30th October, 2011 http://www.timeanddate.com/counters/year2000.html
% Constant parameters
earth_rad=6371000; % in m

%FOR LOOP
v_offset=0;
for sample=1:3
  
v=[0+v_offset 1+v_offset 0+v_offset];
vmag=7200;

[a,inc,RAAN,w,M0,e,R_magnitude,v,r]=orbital_elements(v,vmag);
v=v*(vmag/sqrt(dot(v,v)));
%angular velocity of chief satellite
angular_v=sqrt(mu/norm(R_magnitude)^3);
%sat 2 state vector in meter 
r1=[7695 0 2]*1e3;
%relative initial position vector in LVLH frame

delta_r=r1-r;
delta_r(1)=1000;
%intial condition for bounded relative motion in LVLH frame
vy1=-(delta_r(1))*2*angular_v;
delta_r(3)=-2*delta_r(1);
%VELOCITY vector in LVLH frame for sat 2 
v1=[0 vy1 0];


%initial state in LVLH frame
omega=RAAN;
Inc=inc; 
%direction cosine matrix for  ECI-LVLH coordinates transform
  DCM=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
    -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
    sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
%initial relative state in LVLH frame
initial_state=delta_r;
%calculate parameter for HCW linearized solution
phase_x=atan((angular_v*v1(1))/(3*delta_r(1)+2*v1(2)/angular_v));
phase_z=atan(-(1/angular_v)*(v1(3)/delta_r(3)));
A_z=initial_state(3)/real(cos((phase_z)));
A_x=-(3*initial_state(1)+2*v1(2)/angular_v)/cos(real(phase_x));
%initial relative state in ECI frame
initial_stateECI=inv(DCM)*initial_state';
initial_vstateECI=inv(DCM)*v1';
intial_relativeECI=[initial_stateECI;initial_vstateECI];
%calculate velocity vector of sat 2 in ECI frame
v1=v+initial_vstateECI';
%calculate position vector of sat 2 in ECI frame
r1=r+initial_state;
%velocity magnitude of sat 2
vmag1=norm(v1);
%normalized velocity vector
 v1=v1/norm(v1);
%  [a1,Inc1,RAAN1,w1,M1,e1,R_magnitude1,v1,r1]=orbital_elements(v1,vmag1);
[a1,Inc1,RAAN1,w1,M1,e1,v1,r1]=orbital_elements1(v1*vmag1,r1);
w1=0;
%velocity magnitude of deputy satellite
v1_mag=norm(v1);
v_offset=v_offset;  
%Plotting the ECI Axes
lim=(1-e)*a;
clf
axis([-lim, lim, -lim, lim, -lim, lim])	
view(150,15) 
axis equal
shg
hold on
grid on
% title('Comparison between Predicted Satellite Eliptical Orbit with Extended Kalman Filter and True orbit  ');
line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')
%plot the intial velocity and position vector
%line([r(1)+v(1)*1000],[r(2)+v(2)*1000],[r(3)+v(3)*1000]);
%line(r(1),r(2),r(3));

%Plotting the Earth
equat_rad=6378137.00;
polar_rad=6356752.3142;
% [xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
% load('topo.mat','topo','topomap1');
%     topo2 = [topo(:,181:360) topo(:,1:180)];
%     pro.FaceColor= 'texture';
%     pro.EdgeColor = 'none';
%     pro.FaceLighting = 'phong';
%     pro.Cdata = topo2;
%    earth= surface(xx,yy,zz,pro)
%     colormap(topomap1)
omega_earth = 7.292115855377074e-005; % (rad/sec)  
Go = 1.727564365843028; % (rad)  
GMST = Go + omega_earth*86400*(J2000_days + 0.5);
GMST = GMST - 2*pi*floor(GMST/(2*pi));
GMST_deg=GMST*(180/pi)
% % draw dynamic satellite motion
 orbital_period=norm(sqrt((4*pi*pi*a*a*a)/mu));
 k=0;
 Mean_anomaly=M0;
 for time_instant=1:rev_number*ceil(orbital_period/satrunspeed)
    
    k=k+1;
 %calculate angular velocity of both satellites same for both due to same
 %semi-major axis
%    angular_v=180*orbital_period/satrunspeed/pi; 
   angular_v=sqrt(mu/a^3);
    %calculate eccentric anomaly from mean anomaly
    E_new=real(Mean_anomaly);
    for i=1:5
         E_new=real(Mean_anomaly+(Mean_anomaly + e*sin(E_new) - E_new)/(1 - e*cos(E_new)));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    true_anomaly=2*atan(sqrt((1+e)/(1-e))*tan(E_new/2));
    omega=RAAN;
    Inc=inc; 
%compute radius component r
    R=a*(1-e*sin(E_new));
    R1=a1*(1-e1*sin(E_new));
    rint(k)=a*(1-e*sin(E_new));
    R_magnitude=R;
    Xeci(sample,k) = R_magnitude*(cos(w + true_anomaly)*cos(omega) - sin(w+true_anomaly)*sin(omega)*cos(Inc));
    Yeci(sample,k)= R_magnitude*(cos(w + true_anomaly)*sin(omega) + sin(w+true_anomaly)*cos(omega)*cos(Inc));
    Zeci(sample,k) = R_magnitude*(sin(w + true_anomaly)*sin(Inc));
    %ECI coordinate of sat 2
    RAAN1=0.002;
    Xeci1(sample,k) = R1*(cos(w1 + true_anomaly)*cos(RAAN1) - sin(w1+true_anomaly)*sin(RAAN1)*cos(Inc1));
    Yeci1(sample,k)= R1*(cos(w1 + true_anomaly)*sin(RAAN1) + sin(w1+true_anomaly)*cos(RAAN1)*cos(Inc1));
    Zeci1(sample,k) = R1*(sin(w1 + true_anomaly)*sin(Inc1));
%relative state vector in ECI frame
relative=[Xeci1(sample,k)-Xeci(sample,k) Yeci1(sample,k)-Yeci(sample,k) Zeci1(sample,k)-Zeci(sample,k)];
DCM=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
    -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
    sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
relative_distance=DCM*relative';
    %satelite 1 position vector
sat1_position=[Xeci;Yeci;Zeci];
%% HCW relative distance
x_t=A_x*cos(M0+phase_x)+4*initial_state(1);
y_t=-2*A_x*sin(M0+phase_x)+2*A_x*sin(phase_x);
z_t=A_z*cos(M0+phase_z);
%relative distance from HCW
relative_HCW(k)=norm([x_t y_t z_t]);

   %% Projected Circular Orbit
   %distance parameter of intersatellite network
  distance_parameter=A_x;
  %determine initial phase from initial condition constraint and distance
  %paramter
  intial_phase1=asin(real(2*initial_state(1)/distance_parameter));
%satellite displacement and velcoity components in LVLH frame 
relative_x(k)=-0.5*distance_parameter*cos(real(true_anomaly+intial_phase1));
relative_xvelocity=0.5*distance_parameter*angular_v*sin(real(true_anomaly+intial_phase1));
relative_y(k)=distance_parameter*sin(real(true_anomaly+intial_phase1));
relative_yvelocity=distance_parameter*angular_v*cos(real(true_anomaly+intial_phase1));
relative_z(k)=distance_parameter*cos(real(true_anomaly+intial_phase1));
relative_zvelocity=distance_parameter*angular_v*sin(real(true_anomaly+intial_phase1));
relative_v1=[relative_xvelocity relative_yvelocity relative_zvelocity];

%projected circular orbit solution
   projected_circular=[relative_x(k) relative_y(k) relative_z(k)];
magnitudepro(k)=norm(projected_circular);   
   %% General Circular Orbit
   %distance parameter of intersatellite network
  distance_parameter=A_x;
  %determine initial phase from initial condition constraint and distance
  %paramter
intial_phase1=asin(real(2*r1(1)/distance_parameter));
relative_x1(k)=-0.5*distance_parameter*cos(real(true_anomaly+intial_phase1));
relative_xvelocity1=0.5*distance_parameter*angular_v*sin(true_anomaly+intial_phase1);
relative_y1(k)=distance_parameter*sin(real(true_anomaly+intial_phase1));
relative_yvelocity1=distance_parameter*angular_v*cos(real(true_anomaly+intial_phase1));
relative_z1(k)=(sqrt(3)/2)*distance_parameter*cos(real(true_anomaly+intial_phase1));
relative_zvelocity1=(sqrt(3)/2)*distance_parameter*angular_v*sin(true_anomaly+intial_phase1);
relative_v11=[relative_xvelocity relative_yvelocity relative_zvelocity];
%projected circular orbit coordinates
general_circular=[relative_x1(k) relative_y1(k) relative_z1(k)];
%inter-satellite distance
magnitudegen(k)=norm(general_circular);   
   %%
   
    %check whether satellite is dropping down to earth
    if abs(rint(k))<=earth_rad
       error='satellite drops';
       return;
    end
%     if sample>=1
%     M0=true_anomaly;
    DCM=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
    -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
    sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
    %ECI relative coordinates
%     ECI_relative=[Xeci(sample,k)-Xeci(1,k);Yeci(sample,k)-Yeci(1,k);Zeci(sample,k)-Zeci(1,k)];
%     LVLH_coordinate=DCM*ECI_relative;
    %transform projected circular orbit solution into ECI frame
    ECI_projected=inv(DCM)*projected_circular';
    %transform geneal circular orbit solution into ECI frame
    ECI_projected_general=inv(DCM)*general_circular';
    %ECI represetation of deputy satellite
    ECI_sat2=ECI_projected+[Xeci(sample,k) Yeci(sample,k) Zeci(sample,k)]';
    ECI_sat2_x(k)=ECI_sat2(1);
    ECI_sat2_y(k)=ECI_sat2(2);
    ECI_sat2_z(k)=ECI_sat2(3);
    ECI_sat2_general=ECI_projected_general+[Xeci(sample,k) Yeci(sample,k) Zeci(sample,k)]';
    ECI_sat2_generalx(k)=ECI_sat2_general(1);
    ECI_sat2_generaly(k)=ECI_sat2_general(2);
    ECI_sat2_generalz(k)=ECI_sat2_general(3);
    %position vector components in LVLH frame
%     radial(sample,k)=LVLH_coordinate(1);
%     along(sample,k)=LVLH_coordinate(2);
%     cross(sample,k)=LVLH_coordinate(3);
%     %calculate norm of position in LVLH frame
%     lvlh_magnitude(sample,k)=norm(LVLH_coordinate);
%     end
 Mean_anomaly=norm(Mean_anomaly+sqrt(mu/(a*a*a))*satrunspeed); %Updating Mean Anomaly for next iteration

end
end
%% plot genneral circular orbit 
figure(2)
subplot(3,1,1)
plot3 (Xeci(3,1),Yeci(3,1),Zeci(3,1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','c','MarkerSize', 10);
hold on
plot3 (Xeci(3,2:length(ECI_sat2_x)),Yeci(3,2:length(ECI_sat2_x)),Zeci(1,2:length(ECI_sat2_x)),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
hold on
plot3 (ECI_sat2_generalx(1),ECI_sat2_generaly(1),ECI_sat2_generalz(1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','k','MarkerSize', 10);
hold on
plot3 (ECI_sat2_generalx(2:length(ECI_sat2_generalx)),ECI_sat2_generaly(2:length(ECI_sat2_generalx)),ECI_sat2_generalz(2:length(ECI_sat2_generalx)),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
title('Deputy-Chief GCO Formation Flying in ECI frame')
ylabel('y-coordinate(m) ')
zlabel('z-coordinate(m)')
xlabel('x-coordinate(m)')
legend('Sat 1 initial position','Sat 1 Trajectory in ECI frame','Sat 2 initial position','Sat 2 Trajectory in ECI frame')
subplot(3,1,2)
plot3 (relative_y1,relative_z1,relative_x1,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
ylabel('y-coordinate(m) ')
zlabel('z-coordinate(m)')
xlabel('x-coordinate(m)')
legend('The motion of Deputy Satellite  in Chief-oriented LVLH Frame');
title('General Circular Orbit in LVLH frame')
subplot(3,1,3)
x=linspace(0,length(magnitudepro(1,:))-1,length(magnitudepro(1,:)));
x=x*orbital_period/satrunspeed;
plot (x,magnitudepro/1E3,'LineWidth',2);
hold on
plot (x,magnitudegen/1E3,'LineWidth',2);
xlabel('Time:unit(s)')
ylabel('||r_{LVLH}||(km)')
legend('The inter-satellite distance ||r_{LVLH}|| for Projected Circular Orbit(unit:km)','The inter-satellite distance ||r_{LVLH}|| for General Circular Orbit(unit:km)');
title('Inter-satellite distance ||r_{LVLH}|| Comparison(unit:km)')
% figure(8)
% subplot(3,1,1)
% plot(relative_y1,relative_z1,'LineWidth',2)
% legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
% title('Circular Trajectory projected onto y-z plane in LVLH frame')
% subplot(3,1,2)
% plot(relative_x1,relative_y1,'LineWidth',2)
% legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
% title('Circular Trajectory projected onto x-y plane in LVLH frame')
% subplot(3,1,3)
% plot(relative_x1,relative_z1,'LineWidth',2)
% legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
% title(' Relative Motion on x-z plane in LVLH frame')

%plot ECI motion
% subplot(3,1,3)
% plot(relative_y,relative_z,'LineWidth',2)
% legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
% title('Circular Trajectory projected onto y-z plane in LVLH frame')
%% plot projected circular orbit 
figure(3)
subplot(2,1,1)
plot3 (Xeci(3,3),Yeci(3,3),Zeci(3,3),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','c','MarkerSize', 10);
hold on
plot3 (Xeci(3,2:length(ECI_sat2_x)),Yeci(3,2:length(ECI_sat2_x)),Zeci(1,2:length(ECI_sat2_x)),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
hold on
plot3 (ECI_sat2_x(3),ECI_sat2_y(3),ECI_sat2_z(3),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','k','MarkerSize', 10);
hold on
plot3 (ECI_sat2_x(2:length(ECI_sat2_x)),ECI_sat2_y(2:length(ECI_sat2_x)),ECI_sat2_z(2:length(ECI_sat2_x)),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
ylabel('y-coordinate(m) ')
zlabel('z-coordinate(m)')
xlabel('x-coordinate(m)')
title('Deputy-Chief PCO Formation Flying in ECI frame')
% hold on
%  plot3 (Xeci(4,:),Yeci(4,:),Zeci(4,:),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
% hold on
%  plot3 (Xeci(5,:),Yeci(5,:),Zeci(5,:),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
legend('Deputy Satellite initial position','Deputy Sat Trajectory in ECI frame','Chief Satellite initial position','Chief Sat Trajectory in ECI frame')
subplot(2,1,2)
plot3 (relative_y,relative_z,relative_x,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
ylabel('y-coordinate(m)')
zlabel('z-coordinate(m) ')
xlabel('x-coordinate(m) ')
legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
title('Projected Circular Orbit in LVLH frame')
%plot ECI motion
figure(6)
subplot(3,1,1)
plot(relative_y,relative_z,'LineWidth',2)
legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
title('Circular Trajectory projected onto y-z plane in LVLH frame')
subplot(3,1,2)
plot(relative_x,relative_y,'LineWidth',2)
legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
title('Elliptical Trajectory projected onto x-y plane in LVLH frame')
subplot(3,1,3)
plot(relative_x,relative_z,'LineWidth',2)
legend('The motion of Sat 2 in Sat 1-oriented LVLH Frame');
title('Elliptical Trajectory projected onto x-z plane in LVLH frame')
end