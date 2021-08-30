%% THIS function aims to simulate triangle satellite formation flying for in-plane case and shifted case in LEO 
function triangle_formation
%intial unit velocity vector
v1=[0 0 1];
vmag=7600;
%Gravitational Constant
GM=398.6004418e12;
% Generate a circular orbit
[a,Inc,omega,w,initialMean_anomaly,e,R_magnitude,v1,r1]=orbital_elements(v1,vmag);
%angular velocity of chief satellite
angular_v=sqrt(GM/norm(R_magnitude)^3);
%sat 2 state vector in meter
r2=[7300 0.2 0.2]*1e3;
vy1=-(7.1e6-r1(1))*2*angular_v;
v2=[v1(1) v1(1)+vy1 v1(3)];
vmag1=7500;
v2=v2/vmag1;
%revolution number
rev_number=3;
satrunspeed=300;
%initial position vector in LVLH frame
delta_r=[1e3/6 0 1e3];
%initial velocity vector in LVLH frame
delta_v=[0 -2*angular_v*delta_r(1) 0];
%calculate parameter for HCW linearized solution
phase_x=atan(real((1/angular_v)*(-delta_v(1))/(3*delta_r(1))));
phase_z=atan(-(1/angular_v)*(delta_v(3)/delta_r(3)));
A_z=delta_r(3)/real(cos((phase_z)));
A_x=-(3*delta_r(1)+2*delta_r(2)/angular_v)/cos(real(phase_x));
initial_sat1=[delta_r delta_v]';
%intialization for second satellite projected circular orbits
%with different phases
intial_phase1=0;
intial_phase2=2*pi/3;
intial_phase3=4*pi/3;
delta_r1=[0 0 0];
delta_r1(1)=A_x*cos(intial_phase2)/(-3);
delta_r1(3)=cos(intial_phase2)*A_z;
delta_v1=[0 0 0];
delta_v1(1)=-atan(intial_phase2)*angular_v/(3*delta_r1(1));
delta_v1(3)=-atan(intial_phase2)*angular_v*delta_r1(3);
delta_v1(2)=-2*angular_v*delta_r1(1);
initial_sat2=[delta_r1 delta_v1]';
%initialization for satellite 3
delta_r2=[0 0 0];
delta_r2(1)=A_x*cos(intial_phase3)/(-3);
delta_r2(3)=cos(intial_phase3)*A_z;
delta_v2=[0 0 0];
delta_v2(1)=-atan(intial_phase3)*angular_v/(3*delta_r2(1));
delta_v2(3)=-atan(intial_phase3)*angular_v*delta_r2(3);
delta_v2(2)=-2*angular_v*delta_r2(1);
initial_sat3=[delta_r2 delta_v2]';

%Plotting the ECI Axes
lim=real((1-e)*a/1e5);
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
%observation noise variance range
variacne_measurement=1e2*[2.714e-4 2.7147e-3 2.714e-2 ] ;
for noise_sample=1:length(variacne_measurement)
    
    k=1;
    % % draw dynamic satellite motion
    orbital_period=sqrt((4*pi*pi*a*a*a)/GM);
    %  k=0;
    crosstrackflag=1;
    ECI_flag=1;
    Mean_anomaly=initialMean_anomaly;
    M0=initialMean_anomaly;
    DCM1=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
        -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
        sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
    INTIALR=DCM1*[r1-r2]';
    INTIALV=DCM1*[v1-v2]';
    INTIALR=delta_r'+[0.001 0.001 0.001]';
    INTIALV=delta_v'+[0.001 0.001 0.001]';
    %intial predict error covariance matrix and measurement estimation
    ww=1e-6*[1.8e6, 1.8e6, 1.8e6, 18 ,18 , 18];   
    predicterror_cov_linearunperturb=diag(ww);
    measurment_estimation_linearunperturb=rand([1 6])';
    predicterror_cov_linearunperturb1=diag(ww);
    measurment_estimation_linearunperturb1=rand([1 6])';
    predicterror_cov_linearunperturb2=diag(ww);
    measurment_estimation_linearunperturb2=rand([1 6])';
    %intial state for three satellites
    LOCAL_LVLH=[delta_r delta_v]';
    LOCAL_LVLH1=[delta_r1 delta_v1]';
    LOCAL_LVLH2=[delta_r2 delta_v2]';
    for time_instant=1:rev_number*ceil(orbital_period/satrunspeed)
        
        k=k+1;
        %calculate angular velocity of both satellites same for both due to same
        %semi-major axis
        
        angular_v=sqrt(GM/a^3);
        %calculate eccentric anomaly from mean anomaly
        E_new=real(Mean_anomaly);
        for i=1:5
            E_new=Mean_anomaly+(Mean_anomaly + e*sin(E_new) - E_new)/(1 - e*cos(E_new));
            %         E_new=Mean_anaomaly+e*sin(E_new);
        end
       
        
        true_anomaly=2*atan(sqrt((1+e)/(1-e))*tan(E_new/2));
        
        %compute radius component r
        R=a*(1-e*sin(E_new));
        %direction cosine matrix
  M0=true_anomaly;
  DCM=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
    -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
    sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
angle_flag=1;
%% HCW coordiantes for in-plane triangle formation
     if angle_flag==0
        %for sateliite 1
        x_t=A_x*cos(M0+phase_x)+4*delta_r(1);
        y_t=-2*A_x*sin(M0+phase_x);
        z_t=A_z*cos(M0+phase_z);
          %position vector in LVLh frame
        local_relativePOSITIONLVLH=[x_t y_t z_t]';
          %velcoity components in LVLh frame
        local_relativeVELOCITYLVLH=gradient(local_relativePOSITIONLVLH);
         LOCAL_LVLH=[local_relativePOSITIONLVLH;local_relativeVELOCITYLVLH];
        %eci relative position vector
        ECI_RELATIVE=inv(DCM)*local_relativePOSITIONLVLH;
        LVLH_x(k-1)=local_relativePOSITIONLVLH(1);
        LVLH_y(k-1)=local_relativePOSITIONLVLH(2);
        LVLH_z(k-1)=local_relativePOSITIONLVLH(3);
        %unperturbed relative magnitude
        lvlh_magnitude(k-1)=norm(LOCAL_LVLH);
        %unperturbed relative velocity magnitude
        lvlh_magnitudev(k-1)=norm(gradient(local_relativePOSITIONLVLH));
        %for sateliite 2
        x_t1=A_x*cos(M0+intial_phase2)+4*delta_r(1);
        y_t1=-2*A_x*sin(M0+intial_phase2);
        z_t1=A_z*cos(M0+intial_phase2);
          %position vector in LVLh frame
        local_relativePOSITIONLVLH1=[x_t1 y_t1 z_t1]';
          %velcoity components in LVLh frame
        local_relativeVELOCITYLVLH1=gradient(local_relativePOSITIONLVLH1);
        LOCAL_LVLH1=[local_relativePOSITIONLVLH1;local_relativeVELOCITYLVLH1];
        %eci relative position vector
        ECI_RELATIVE1=inv(DCM)*local_relativePOSITIONLVLH1;
        LVLH_x1(k-1)=local_relativePOSITIONLVLH1(1);
        LVLH_y1(k-1)=local_relativePOSITIONLVLH1(2);
        LVLH_z1(k-1)=local_relativePOSITIONLVLH1(3);
        %unperturbed relative magnitude
        lvlh_magnitude1(k-1)=norm(LOCAL_LVLH1);
        %unperturbed relative velocity magnitude
        lvlh_magnitudev1(k-1)=norm(gradient(local_relativePOSITIONLVLH1));
         %for sateliite 3
        x_t2=A_x*cos(M0+intial_phase3)+4*delta_r(1);
        y_t2=-2*A_x*sin(M0+intial_phase3);
        z_t2=A_z*cos(M0+intial_phase3);
          %position vector in LVLh frame
        local_relativePOSITIONLVLH2=[x_t2 y_t2 z_t2]';
          %velcoity components in LVLh frame
        local_relativeVELOCITYLVLH2=gradient(local_relativePOSITIONLVLH2);
        LOCAL_LVLH2=[local_relativePOSITIONLVLH2;local_relativeVELOCITYLVLH2];
        %eci relative position vector
        ECI_RELATIVE2=inv(DCM)*local_relativePOSITIONLVLH2;
        LVLH_x2(k-1)=local_relativePOSITIONLVLH2(1);
        LVLH_y2(k-1)=local_relativePOSITIONLVLH2(2);
        LVLH_z2(k-1)=local_relativePOSITIONLVLH2(3);
        %unperturbed relative magnitude
        lvlh_magnitude2(k-1)=norm(LOCAL_LVLH2);
        %unperturbed relative velocity magnitude
        lvlh_magnitudev2(k-1)=norm(gradient(local_relativePOSITIONLVLH2));
     end
%% HCW coordiantes for shifted triangle formation
  if angle_flag==1
  %for sateliite 1
        x_t=A_x*cos(M0+phase_x);
        y_t=-2*A_x*sin(M0+phase_x)+2*A_x*sin(intial_phase1);
        z_t=A_z*cos(M0+phase_z); 
        local_relativePOSITIONLVLH=[x_t y_t z_t]';
        %velcoity components in LVLh frame
        local_relativeVELOCITYLVLH=gradient(local_relativePOSITIONLVLH);
         LOCAL_LVLH=[local_relativePOSITIONLVLH;local_relativeVELOCITYLVLH];
        %eci relative position vector
        ECI_RELATIVE=inv(DCM)*local_relativePOSITIONLVLH;
        %LVLH coordiantes for sat 1
        LVLH_x(k-1)=local_relativePOSITIONLVLH(1);
        LVLH_y(k-1)=local_relativePOSITIONLVLH(2);
        LVLH_z(k-1)=local_relativePOSITIONLVLH(3);
        %unperturbed relative magnitude
        lvlh_magnitude(k-1)=norm(local_relativePOSITIONLVLH);
        %unperturbed relative velocity magnitude
        lvlh_magnitudev(k-1)=norm(gradient(local_relativePOSITIONLVLH));
        %LVLH coordiantes for sateliite 2
        x_t1=A_x*cos(M0+intial_phase2);
        y_t1=-2*A_x*sin(M0+intial_phase2)+2*A_x*sin(intial_phase2);
        z_t1=A_z*cos(M0+intial_phase2);
        
        %LVLH POSITION vector for SAT 2 
        local_relativePOSITIONLVLH1=[x_t1 y_t1 z_t1]';
        %LVLH VELOCITY vector for SAT 2 
        local_relativeVELOCITYLVLH1=gradient(local_relativePOSITIONLVLH1);
        LOCAL_LVLH1=[local_relativePOSITIONLVLH1;local_relativeVELOCITYLVLH1];
        %eci relative position vector
        ECI_RELATIVE1=inv(DCM)*local_relativePOSITIONLVLH1;
              %LVLH coordiantes for sat 2
        LVLH_x1(k-1)=local_relativePOSITIONLVLH1(1);
        LVLH_y1(k-1)=local_relativePOSITIONLVLH1(2);
        LVLH_z1(k-1)=local_relativePOSITIONLVLH1(3);
        %unperturbed relative magnitude
        lvlh_magnitude1(k-1)=norm(local_relativePOSITIONLVLH1);
        %unperturbed relative velocity magnitude
        lvlh_magnitudev1(k-1)=norm(gradient(local_relativePOSITIONLVLH1));
        %for sateliite 3
        x_t2=A_x*cos(M0+intial_phase3);
        y_t2=-2*A_x*sin(M0+intial_phase3)+2*A_x*sin(intial_phase3);
        z_t2=A_z*cos(M0+intial_phase3);
        
        local_relativePOSITIONLVLH2=[x_t2 y_t2 z_t2]';
        local_relativeVELOCITYLVLH2=gradient(local_relativePOSITIONLVLH2);
        LOCAL_LVLH2=[local_relativePOSITIONLVLH2;local_relativeVELOCITYLVLH2];
        %eci relative position vector
        ECI_RELATIVE2=inv(DCM)*local_relativePOSITIONLVLH2;
        %LVLH coordiantes for sat 3
        LVLH_x2(k-1)=local_relativePOSITIONLVLH2(1);
        LVLH_y2(k-1)=local_relativePOSITIONLVLH2(2);
        LVLH_z2(k-1)=local_relativePOSITIONLVLH2(3);
        %unperturbed relative magnitude
        lvlh_magnitude2(k-1)=norm(local_relativePOSITIONLVLH2);
        %unperturbed relative velocity magnitude
        lvlh_magnitudev2(k-1)=norm(gradient(local_relativePOSITIONLVLH2));
  end
   %%  plot triangle formation flying in LVLH frame
   if k>=2
   perturbsatellite1=plot3 (LVLH_x(k-1),LVLH_y(k-1),LVLH_z(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
perturbsatellite2=plot3 (LVLH_x1(k-1),LVLH_y1(k-1),LVLH_z1(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
perturbsatellite3=plot3 (LVLH_x2(k-1),LVLH_y2(k-1),LVLH_z2(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 6);
%% commented part plot the line connects satellite position at two consecutive time instances
% line([relative_xperturb relative_x1perturb],[relative_yperturb relative_y1perturb],[relative_zperturb relative_z1perturb],'Color', 'blue', 'LineWidth', 2);
%  line([relative_x1perturb relative_x2perturb],[relative_y1perturb relative_y2perturb],[relative_z1perturb relative_z2perturb],'Color', 'blue', 'LineWidth', 2);
%  line([relative_xperturb relative_x2perturb],[relative_yperturb relative_y2perturb],[relative_zperturb relative_z2perturb],'Color', 'blue', 'LineWidth', 2);
%3D plot for trangle formation in LVLH frame
perturbsatellite11=plot3 (LVLH_x(1),LVLH_y(1),LVLH_z(1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 11);
perturbsatellite22=plot3 (LVLH_x1(1),LVLH_y1(1),LVLH_z1(1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 11);
perturbsatellite33=plot3 (LVLH_x2(1),LVLH_y2(1),LVLH_z2(1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 11);
legend([perturbsatellite1, perturbsatellite2,perturbsatellite3],'Satellite 1','Satellite 2','Satellite 3');
legend([perturbsatellite11, perturbsatellite22,perturbsatellite33,perturbsatellite1, perturbsatellite2,perturbsatellite3],'Satellite 1 Initial Position','Satellite 2 Initial Position','Satellite 3 Initial Position','Satellite 1','Satellite 2','Satellite 3');
xlabel('x-coordinate')
ylabel('y-coordinate')
zlabel('z-coordinate')
   
   end

  %unperturbed relative magnitude
        lvlh_magnitude(k-1)=norm( LOCAL_LVLH(1:3));
        %unperturbed relative velocity magnitude
        lvlh_magnitudev(k-1)=norm( LOCAL_LVLH(4:6));
         %unperturbed relative magnitude
        lvlh_magnitude1(k-1)=norm( LOCAL_LVLH1(1:3));
        %unperturbed relative velocity magnitude
        lvlh_magnitudev1(k-1)=norm( LOCAL_LVLH1(4:6));
         %unperturbed relative magnitude
        lvlh_magnitude2(k-1)=norm( LOCAL_LVLH2(1:3));
        %unperturbed relative velocity magnitude
        lvlh_magnitudev2(k-1)=norm( LOCAL_LVLH2(4:6));

%update mean anomaly
Mean_anomaly=Mean_anomaly+sqrt((GM/(a*a*a)))*satrunspeed;
end
end

%% plot intersatellite distance and relative velocity between sat 1-2 1-3 2-3
figure(7)
subplot(2,1,1)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);

plot(x,lvlh_magnitude,'LineWidth',2);
hold on
plot(x,lvlh_magnitude1,'LineWidth',2);
hold on
plot(x,lvlh_magnitude2,'LineWidth',2);
grid on
xlabel('Time(unit:s)');
legend('||r_{true}|| for Sat 1','||r_{true}|| for Sat 2','||r_{true}|| for Sat 3')
subplot(2,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,lvlh_magnitudev,'LineWidth',2);
hold on
plot(x,lvlh_magnitudev1,'LineWidth',2);
hold on
plot(x,lvlh_magnitudev2,'LineWidth',2);
grid on
xlabel('Time(unit:s)');
legend('Velocity Magnitude ||v_{Sat1}|| in LVLH frame','Velocity Magnitude ||v_{Sat2}|| in LVLH frame','Velocity Magnitude ||v_{Sat3}|| in LVLH frame')

end