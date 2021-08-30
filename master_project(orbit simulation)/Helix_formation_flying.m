%% This program aims to apply Kalman filter for deputy satellite motion 
%% state estimation in LEO under different cases, evaluate its tracking performance 
%% and visualize deputy satellite motion with/without perturbation
function Helix_formation_flying
%intial velocity vector
v1=[0 0 1];
vmag=7889;
%Gravitational Constant
GM=398.6004418e12;
%generate a circular chief satellite orbit
[a,Inc,omega,w,initialMean_anomaly,e,R_magnitude,v1,r1]=orbital_elements(v1,vmag);
%angular velocity of chief satellite
angular_v=sqrt(GM/norm(R_magnitude)^3);
%deputy satellite state vector in meter
r2=[7300 0.2 0.2]*1e3;
vy1=-(7.1e6-r1(1))*2*angular_v;
v2=[v1(1) v1(1)+vy1 v1(3)];
vmag1=7500;
v2=v2/vmag1;
%initial position vector in LVLH frame
delta_r=[1e3/6 0 1e3];
%initial velocity vector in LVLH frame
delta_v=[0 -2*angular_v*delta_r(1) 0];
%calculate parameter for HCW linearized solution
phase_x=atan(real((1/angular_v)*(-delta_v(1))/(3*delta_r(1))));
phase_z=atan(-(1/angular_v)*(delta_v(3)/delta_r(3)));
A_z=delta_r(3)/real(cos((phase_z)));
A_x=-(3*delta_r(1)+2*delta_r(2)/angular_v)/cos(real(phase_x));
[a1,Inc1,omega1,w1,Mean_anomaly1,e1,R_magnitude1,v2,r2]=orbital_elements(v2,vmag1);
%Number of Days since Jan 1, 2000
J2000_days=103752/24;
%satellite speed
satrunspeed=200; 
%revolution time
rev_number=3;
%earth radius
EARTH_RAIDUS=6371000;% unit:m
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
%% Plotting the Earth this commented part is not mandatory for plotting, depends on user demand
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
GMST_deg=GMST*(180/pi);%2.714e-5 2.714e-4 2.714e-3 2.714e-2 2.714e-1 2.714e0
%orbital period
orbital_period=sqrt((4*pi*pi*a*a*a)/GM);
%% averaging for different simulation samples
%observation noise variance array
 variacne_measurement=1e8*[2.714e-6 2.7147e-5 2.714e-4 ];
%total number of simulations
simul_number=1;
%intial array to store estimation error array for unperturbed with J2 and
%with air drag forces 
final_total_errordragunperturb=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
final_total_errordragunperturb1=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
final_total_errordrag=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
final_total_errordrag1=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
final_total_errorj2=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
final_total_errorj21=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
% initial averaged residual array to store residual norm values
finalredisualdrag_norm=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finalredisualdrag_normvelocity=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finalredisualunper_normr=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finalredisualunper_normv=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finalredisualj2_normr=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finalredisualj2_normv=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finallvlh_magnitude=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finallvlh_magnitudev=zeros(1,rev_number*ceil(orbital_period/satrunspeed));
finalunperturbpos=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
finalunperturbvel=zeros(length(variacne_measurement),rev_number*ceil(orbital_period/satrunspeed));
%loop over all simulation samples
for local_simu_number=1:simul_number
%loop over all observation noise variances
for noise_sample=1:length(variacne_measurement)
    k=1;
 M0=initialMean_anomaly;
 DCM1=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
    -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
    sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
 %INITIAL state vector in ECI frame
 initial_state=[inv(DCM1)*delta_r';inv(DCM1)*delta_v']+[r1';v1'];
 INTIALR=delta_r'+[0.001 0.001 0.001]';
 INTIALV=delta_v'+[0.001 0.001 0.001]';
%intial predict error covariance matrix and measurement estimation for
%classic HCW model and drag model
ww=1e-6*[1.8e6, 1.8e6, 1.8e6, 18 ,18 , 18];
predicterror_cov_linear=1e-6*diag(ww);
measurment_estimation_linear=rand([1 6])';
predicterror_cov_linearunperturb=diag(ww);
measurment_estimation_linearunperturb=rand([1 6])';
predicterror_cov_linearunperturb1=diag(ww);
measurment_estimation_linearunperturb1=[INTIALR;INTIALV];
%initial state under J2
J2_systemvalue=[delta_r';delta_v'];
%initial state guess and predict error covariance matrix for J2 model
measurment_estimation_J2=rand([1 6])';
predicterror_cov_linearJ2=diag(ww);
Mean_anomaly=0;
%initial relative distance
relative_rmag(1)=norm(measurment_estimation_linear(1:3));
relative_magnitude=[INTIALR;INTIALV];
%intial drag state
relative_state_drag=[delta_r';delta_v'];
LOCAL_J2PVPH=[delta_r';delta_v'];
%loop over different orbitial periods
 for time_instant=1:rev_number*ceil(orbital_period/satrunspeed)
    k=k+1;
 %calculate angular velocity of both satellites same for both due to same
 %semi-major axis
     angular_v=sqrt(GM/a^3);
    %calculate eccentric anomaly from mean anomaly
    E_new=real(Mean_anomaly);
    for i=1:5
         E_new=Mean_anomaly+(Mean_anomaly + e*sin(E_new) - E_new)/(1 - e*cos(E_new));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    E_new1=Mean_anomaly1;
    for i=1:5
         E_new1=Mean_anomaly1+(Mean_anomaly1 + e*sin(E_new1) - E_new1)/(1 - e*cos(E_new1));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
   %calculate true anomaly from calculated eccentric anomaly
    true_anomaly=2*atan(sqrt((1+e)/(1-e))*tan(E_new/2));
    true_anomaly1=2*atan(sqrt((1+e)/(1-e))*tan(E_new1/2));
    %compute radius component r
    R=a*(1-e*sin(E_new));
    rint(k)=a*(1-e*sin(E_new));
    R1=a1*(1-e1*sin(E_new1));
    r11(k)=a*(1-e1*sin(E_new1));
   %check whether satellite is dropping down to earth
    if rint(k)<=EARTH_RAIDUS
       error='satellite drops'
       return;
   end
   if r11(k)<=EARTH_RAIDUS
       error='satellite drops'
       return;
   end
   
   % ECI coordinates presentation without perturbation
    Xeci = R*(cos(w + true_anomaly)*cos(omega) - sin(w+true_anomaly)*sin(omega)*cos(Inc));
    Yeci= R*(cos(w + true_anomaly)*sin(omega) + sin(w+true_anomaly)*cos(omega)*cos(Inc));
    Zeci = R*(sin(w + true_anomaly)*sin(Inc));
%calculate orbital height of chief satellite
orbital_heightchief(k-1)=norm([Xeci Yeci Zeci])-EARTH_RAIDUS;
    %ECI coordinates of deputy satellite
    Xeci1 = R1*(cos(w1 + true_anomaly)*cos(omega1) - sin(w1+true_anomaly)*sin(omega1)*cos(Inc1));
    Yeci1= R1*(cos(w1 + true_anomaly)*sin(omega1) + sin(w1+true_anomaly)*cos(omega1)*cos(Inc1));
    Zeci1 = R1*(sin(w1 + true_anomaly)*sin(Inc1));
    xcor1(k-1)=Xeci1;
    ycor1(k-1)=Yeci1;
    zcor1(k-1)=Zeci1;
relative_r=[Xeci1-Xeci;Yeci1-Yeci;Zeci1-Zeci];
relative_rmag(k)=norm(relative_r);

%% Linearized solution for state transition matrix for relative position and velocity of deputy satellite
%% with respect to chief 
Position_transferXX=[(4-3*cos(real(Mean_anomaly))) 0 0;
     (6*sin(real(Mean_anomaly))-6*real(Mean_anomaly)) 1 0;
     0 0 cos(real(Mean_anomaly))];
 
 Velocity_transferXV=[(sin(real(Mean_anomaly))/angular_v)  2*(1-cos(real(Mean_anomaly)))/angular_v  0;
     2*(cos(real(Mean_anomaly))-1)/angular_v  (4*sin(real(Mean_anomaly))-3*real(Mean_anomaly))/angular_v  0;
     0 0 sin(real(Mean_anomaly))/angular_v];
 Position_transferVX=[3*angular_v*sin(Mean_anomaly) 0 0;
     -6*angular_v*(1-cos(Mean_anomaly))  0 0;
     0 0 -sin(angular_v)*angular_v];
 Velocity_transferVV=[cos(real(Mean_anomaly)) 2*sin(real(Mean_anomaly)) 0;
     -2*sin(real(Mean_anomaly)) -3+4*cos(real(Mean_anomaly)) 0;
     0 0 cos(real(Mean_anomaly))];
 %state transition matrix under linear solution
 Newstate=[Position_transferXX,Velocity_transferXV;
    Position_transferVX,Velocity_transferVV];

  
%% VELOCITY vector for both satellites in ECI frame
    VX1 = angular_v*(-sin(w1 + true_anomaly1)*cos(omega1) - angular_v*cos(w1+true_anomaly1)*sin(omega1)*cos(Inc));
    VY1= angular_v*(-sin(w1 + true_anomaly1)*sin(omega1) + angular_v*cos(w1+true_anomaly1)*cos(omega1)*cos(Inc));
    VZ1 = angular_v*(cos(w1 + true_anomaly1)*sin(Inc));
    v_vector1=R*[ VX1; VY1; VZ1];
    VX2 = -sin(w + true_anomaly)*cos(omega) - cos(w+true_anomaly)*sin(omega)*cos(Inc);
    VY2= (-sin(w + true_anomaly)*sin(omega) + cos(w+true_anomaly)*cos(omega)*cos(Inc));
    VZ2 = (cos(w + true_anomaly)*sin(Inc));
    v_vector2=R1*[ VX2; VY2; VZ2];
    %true state
    true_state=[Xeci1;Yeci1;Zeci1;v_vector2];   
    %ECI-relative velocity vector
    relative_v=[v_vector2-v_vector1];
    
%% CALCULATE projection polar coordinates on the earth surface
%calculate latitude
    lat(k)=atan(Zeci/sqrt(Xeci*Xeci+Yeci*Yeci));    
    lat1(k)=atan(Zeci1/sqrt(Xeci1*Xeci1+Yeci1*Yeci1))*(180/pi);
    %calculate longitude
     ECIX=[cos(GMST) sin(GMST) 0];
    Pos=[real(Xeci) real(Yeci) 0];
    cvec = cross(ECIX,Pos);
    Pos1=[real(Xeci1) real(Yeci1) 0];
    cvec1 = cross(ECIX,Pos1);
    angleyz = mod(sign(dot([0 0 1],cvec))*atan2(norm(cvec),dot(ECIX,Pos)),2*pi);
    long(k) =(180/pi)* angleyz;    
    angleyz1 = mod(sign(dot([0 0 1],cvec1))*atan2(norm(cvec1),dot(ECIX,Pos1)),2*pi);
    long1(k) =(180/pi)* angleyz1;
 
   
 %%
%relative state combine perturbation
     relative_state=[relative_r;relative_v];
%direction cosine matrix
  M0=true_anomaly;
  DCM=[cos(omega)*cos(w+M0)-sin(omega)*sin(w+M0)*cos(Inc) sin(omega)*cos(w+M0)+cos(omega)*sin(w+M0)*cos(Inc) sin(w+M0)*sin(Inc);
    -cos(omega)*sin(w+M0)-sin(omega)*cos(w+M0)*cos(Inc) -sin(omega)*sin(w+M0)+cos(omega)*cos(w+M0)*cos(Inc) cos(w+M0)*sin(Inc);
    sin(omega)*sin(Inc) -cos(omega)*sin(Inc) cos(Inc)];
  
  %% HCW relative distance
x_t=A_x*cos(M0+phase_x)+4*delta_r(1);
y_t=-2*A_x*sin(M0+phase_x)+2*A_x*sin(phase_x);
z_t=A_z*cos(M0+phase_z);
vx_t=-A_x*angular_v*sin(M0+phase_x);
vy_t=-2*A_x*angular_v*cos(M0+phase_x);
vz_t=-A_z*angular_v*sin(M0+phase_z); 
local_relativePOSITIONLVLH=[x_t y_t z_t]';
local_relativeVELOCITYLVLH=gradient(local_relativePOSITIONLVLH);
lvlh_magnitudev(k-1)=norm(local_relativeVELOCITYLVLH);
%eci relative position vector
ECI_RELATIVE=inv(DCM)*local_relativePOSITIONLVLH;
%DEPUTY satellite coordinates in LVLH frame
LVLH_x(k-1)=local_relativePOSITIONLVLH(1);
LVLH_y(k-1)=local_relativePOSITIONLVLH(2);
LVLH_z(k-1)=local_relativePOSITIONLVLH(3);  
%unperturbed relative magnitude
lvlh_magnitude(k-1)=norm(local_relativePOSITIONLVLH);
%relative state vector of deputy satellite in LVLH frame
LOCAL_LVLH=[local_relativePOSITIONLVLH;local_relativeVELOCITYLVLH];

%% kalman filter implementation
J2perturbation_matrix=[0 0 0 0 0 0]';
atmospheric_drag_term=[0 0 0 0 0 0]';
%create zero matrix
zero_matrix=zeros(3);
%creat identity matrix
Identity=eye(3);
% k matrix 
k_matrix=[angular_v^2 0 0;
    0 angular_v^2 0;
    0 0 0];
omeag_matrix=[0 2*angular_v 0;
    -2*angular_v 0 0;
    0  0 0];
%establish j matrix
v=[2*angular_v^2, -angular_v^2,-angular_v^2];
J=diag(v);
%process noise covariance matrix
measurement_noiselinear=1e9*diag([2.7e-7;2.7e-7;2.7e-7;2.7e-7;2.7e-7;2.7e-7]);
%random noise to generate true state
process_noiselinear=sqrt(2.7e2)*rand(1,6)';
%measurement noise covariance 
measurement_noise_linear=sqrt(variacne_measurement(noise_sample))*randn(1,6)';
%establish system matrix 
 systemmatrixunperturb=([zero_matrix Identity;
     k_matrix+J  omeag_matrix ])+eye(6);
% true state with random noise
relative_magnitude=systemmatrixunperturb^(k-1)*[delta_r';delta_v']+process_noiselinear;
  model_state=relative_magnitude;
%position and velocity magnitude 
lvlh_magnitudev(k-1)=norm(model_state(4:6));
lvlh_magnitude(k-1)=norm(model_state(1:3));
%ECI Coordinate represenation of deputy saetllite
eci_vector=norm([Xeci Yeci Zeci]')+lvlh_magnitude(k-1);
orbital_heightdeputy(k-1)=norm(eci_vector)-EARTH_RAIDUS;
    % measurement_noiselinear=diag([distanceerror;velocityerror]);
%% modified state transition model with atmospheric drag force model
%h:integral constant r;raidus of deputy satellite u:gravitational constant
h=2.12e15;
u=398.6004418e12;
drag_constant=2.2e-3;
%a flag to control if we wanna to aplly state estimation for drag model
drag_flag=1;
if drag_flag==1
    lower_left=[0 0 0; 0 3+12*drag_constant^2 0; 0 0 -1];
    lower_right=[0 2 0; -2 0 0; 0 0 0];
Atmospheric_HCW=real([zeros(3) eye(3);
    lower_left lower_right]);
%due to nonlinear state transition model, jacobian matrix is calculated (3*u*h^(-2)*R_magnitude^3)/exp(-2*drag_constant*lat(k))
%
systemmatrix=Atmospheric_HCW+eye(6);
 %establish linear observation model
linear_obervationmodel=eye(6);
%process noise for linear model 1e7*[2.714e-7 3.7147e-7 4.714e-7 ]
variacne_measurementdrag=1e9*[2.714e-5 2.7147e-4 2.714e-3 ];
R_matrix_lineardrag=variacne_measurementdrag(noise_sample)*eye(6);
%measurement noise matrix
measurement_noiselineardrag=1e7*diag([2.7e-7;2.7e-7;2.7e-7;2.7e-7;2.7e-7;2.7e-7]);
% measurement_noiselinear=distance_error(k-1);
measurement_noise_lineardrag=sqrt(variacne_measurementdrag(noise_sample))*randn(1,6)';
R_matrix_linear=variacne_measurement(noise_sample)*eye(6);
R_matrix_linearun=variacne_measurement(noise_sample)*eye(6);

%% construct HCW model under atmospheric drag force
if drag_flag==1
    %true anomaly for initial circular orbits
    time=true_anomaly;
    initial_state=[r1-r2,v1-v2]';
    %drag constant
    kepa=sqrt(abs(1-12*drag_constant^2));
    %system matrix formulation
   upperleft=[1 4*sin(kepa*time)/kepa+(1-4/kepa^2)*time  2*(1-4/kepa^2)*(sin(kepa*time)/kepa-time);
       0 4*(cos(kepa*time)-1)/kepa^2+1 2*(1-4/kepa^2)*(cos(kepa)-1);
       0 2*(cos(kepa*time)-1)/kepa^2 (1-4/kepa^2)*cos(kepa*time)+4/kepa^2];
   upperright=[2*(1-cos(kepa*time))/kepa^2 0 0;
       2*sin(kepa*time)/kepa 0 0;
       sin(kepa*time)/kepa 0 0];
   lowerleft=[0 -2*sin(kepa*time)/kepa -kepa*(1-4/kepa^2)*sin(kepa*time);
       0 0 0 ;
       0 0 0];
   lowerright=[cos(kepa*time) 0 0;
       0 cos(time) sin(time)
       0 -sin(time) cos(time)];
   final_statetransitionmodel=[upperleft upperright
       lowerleft lowerright];
   %satellite state under drag force  
    process_noiselineardrag=awgn(randn(1,6)',1e9*[2.7e-7]);
   FINAL_DRAG=final_statetransitionmodel*[delta_r';delta_v']+process_noiselineardrag;
   %drag-modified coordinates in LVLH frame
   drag_x(k-1)=relative_state_drag(1);
   drag_y(k-1)=relative_state_drag(2);
   drag_z(k-1)=relative_state_drag(3);
  
%% Drag model state
%true state with random noise
process_noiselineardrag=awgn(randn(1,6)',1e9*[2.7e-7]);
 relative_state_drag=systemmatrix^(k-1)*[delta_r';delta_v']+process_noiselineardrag;
%position vector magnitude 
drag_magnitude(k-1)=norm(FINAL_DRAG(1:3));
%velocity vector magnitude
drag_magnitudev(k-1)=norm(FINAL_DRAG(4:6));
if drag_flag==1
y_linear=linear_obervationmodel*relative_state_drag+measurement_noise_lineardrag;
end 
%position and velocity estimation error under drag model trace values of Predict error covariance
if k==2
error=norm(measurment_estimation_linear(1:3)-relative_state_drag(1:3));
total_errordrag(noise_sample,1)=real(20*log(error/norm(FINAL_DRAG(1:3))));
 errordrag=norm(measurment_estimation_linear(4:6)-relative_state_drag(4:6));
total_errordrag1(noise_sample,1)=real(20*log(errordrag/norm(FINAL_DRAG(4:6))));
dragtrace(noise_sample,1)=trace(predicterror_cov_linear(1:3,1:3));
dragtrace1(noise_sample,1)=trace(predicterror_cov_linear(4:6,4:6));
   end
%predicted error covariance matrix estimate for linear model with drag
%force
predict_error_cov=systemmatrix*predicterror_cov_linear*systemmatrix'+measurement_noiselineardrag;
%predicted state estimate
estimate_measurment=systemmatrix*measurment_estimation_linear;
%calculate residual based on linear observation model
residual_linear=y_linear-linear_obervationmodel*(estimate_measurment+J2perturbation_matrix+atmospheric_drag_term);
% residual norm for with drag force
redisualdrag_norm(noise_sample,k-1)=norm(residual_linear(1:3));
%norm of velocity portion of residual
redisualdrag_normvelocity(noise_sample,k-1)=norm(residual_linear(4:6));
%intermediate step to calculate kalman gain
b2=(linear_obervationmodel*predict_error_cov*transpose(linear_obervationmodel)+R_matrix_lineardrag);
%calcualte kalman gain
kalman_gain_linear= predict_error_cov*transpose(linear_obervationmodel)*pinv(b2);
%calculate state estimation with J2 perturbation and atmospheric drag force for linear model 
measurment_estimation_linear=estimate_measurment+kalman_gain_linear*residual_linear;
a_postestimatedrag=measurment_estimation_linear;

%difference between a prior state estimate and true state
drag_residualvar=(relative_state_drag-a_postestimatedrag);
% d_metric represenation for drag model
drag_distance(noise_sample,k-1)=(norm(drag_residualvar(1:3)));
drag_velocity(noise_sample,k-1)=norm(drag_residualvar(4:6));
%calculate updated predict error covariance matrix update for linear model with
%drag force
predicterror_cov_linear=(eye(6)-kalman_gain_linear*linear_obervationmodel)*predict_error_cov;
%trace for position and velocity for drag model
if k>2
%trace for position part
dragtrace(noise_sample,k-1)=trace(predicterror_cov_linear(1:3,1:3));
%trace for velocity
dragtrace1(noise_sample,k-1)=trace(predicterror_cov_linear(4:6,4:6));
end
end
%relative distance with perturbation of drag force
perturbed_intersatellite_distance(k-1)=norm(relative_state_drag(1:3));
end
%% kalman filter evaluation for linear unpuerturbed HCW model
% calculate initial position and velocity estimation error and trace 
 if k==2
error1=measurment_estimation_linearunperturb(1:3)- model_state(1:3);
linear=norm(error1)/norm(LOCAL_LVLH(1:3));
total_errordragunperturb(noise_sample,1)=real(20*log(norm(error1)/norm(model_state(1:3))));
sumdiagonal(noise_sample,1)=trace(predicterror_cov_linearunperturb(1:3,1:3));
traceval(noise_sample,1)=trace(predicterror_cov_linearunperturb(1:3,1:3));
error11=measurment_estimation_linearunperturb(4:6)- model_state(4:6);
total_errordragunperturb1(noise_sample,1)=real(20*log(norm(error11)/norm(model_state(4:6))));
sumdiagonal1(noise_sample,1)=trace(predicterror_cov_linearunperturb(4:6,4:6));
trace1(noise_sample,1)=trace(predicterror_cov_linearunperturb(4:6,4:6));
 end
%calculate measurement value
 y_linearunperturb=linear_obervationmodel* model_state+measurement_noise_linear;
 %predicted error covariance matrix estimate for linear model
 predict_error_covunperturb=systemmatrixunperturb*predicterror_cov_linearunperturb*systemmatrixunperturb'+measurement_noiselinear;
 estimate_measurmentunperturb=systemmatrixunperturb*measurment_estimation_linearunperturb;
%calculate residual based on linear observation model
 residual_linearunperturb=y_linearunperturb-linear_obervationmodel*(estimate_measurmentunperturb);
%norm of residual
redisual_normposition(k-1)=norm(residual_linearunperturb(1:3));
redisual_normvelocity(k-1)=norm(residual_linearunperturb(4:6));
%calculate kalman gain
A=(linear_obervationmodel*predict_error_covunperturb*transpose(linear_obervationmodel)+R_matrix_linear);
 %calcualte kalman gain
kalman_gain_linearunperturb= predict_error_covunperturb*transpose(linear_obervationmodel)*pinv(A);
%a posterior state estimate
a_postestimate=measurment_estimation_linearunperturb;
%calculate state estimation with J2 perturbation and atmospheric drag force for linear model 
measurment_estimation_linearunperturb=estimate_measurmentunperturb+J2perturbation_matrix+atmospheric_drag_term+kalman_gain_linearunperturb*residual_linearunperturb;
%difference between a prior state estimate and true state
%a posterior estimation error
unperturb_residualvar=systemmatrixunperturb*(model_state-a_postestimate)+process_noiselinear;
unperturbvar=(model_state-a_postestimate);
% unperturb_residualvar=measurement_noise_linear+process_noiselinear;
unperturbpos(noise_sample,k-1)=norm(unperturb_residualvar(1:3));
unperturbvel(noise_sample,k-1)=norm(unperturb_residualvar(4:6));
Aposterior_pos(noise_sample,k-1)=norm(unperturbvar(1:3));
Aposterior_vel(noise_sample,k-1)=norm(unperturbvar);
%calculate predict error covariance matrix update for linear model
predicterror_cov_linearunperturb=(eye(6)-kalman_gain_linearunperturb*linear_obervationmodel)*predict_error_covunperturb;
%TRACE of predict error covariance matrix
if k>2
sumdiagonal(noise_sample,k-1)=trace(predicterror_cov_linearunperturb(1:3,1:3));
%trace value for unperturb model
traceval(noise_sample,k-1)=trace(predicterror_cov_linearunperturb(1:3,1:3));
sumdiagonal1(noise_sample,k-1)=trace(predicterror_cov_linearunperturb(4:6,4:6));
%trace value for unperturb model
trace1(noise_sample,k-1)=trace(predicterror_cov_linearunperturb(4:6,4:6));
end

%% J2 perturbed extended kalman filter
 %J2 acceleration components in LVLH frame
 J2=1.083e-3;
 j2_factor=3*GM*J2*EARTH_RAIDUS^2/(2*R^4);
 %angle theta
 theta=w+true_anomaly;
 j2_acceleration=j2_factor*[0;0;0;-1+3*sin(Inc)^2*sin(theta)^2;-2*sin(Inc)^2*sin(theta)*cos(theta);-2*sin(Inc)*cos(theta)*sin(theta)];
 %state vector in LVLH frame under J2 solved by J2 modified HCW
 %equation
 J2_term=3*J2*EARTH_RAIDUS^2/(2*R);
 J2_radial(k-1)=J2_term*((1/3)*sin(Inc)^2*cos(theta)^2+((1/3)*sin(Inc)^2-1)+(1-(2/3)*sin(Inc)^2)*cos(theta));
 J2_ALONG(k-1)=J2_term*((2/3)*sin(Inc)^2*cos(theta)+(1-(2/3)*sin(Inc)^2)*cos(theta))+(j2_factor*sin(Inc)^2*cos(2*theta))/(4*angular_v^2);
J2_cross(k-1)=j2_factor*sin(Inc)*sin(2*theta)/(-5*angular_v^2);
%position vector in LVLH frame
 relative_J2LVLH=[J2_radial(k-1);J2_ALONG(k-1);J2_cross(k-1)];
 %modified J2 perturbed HCW model
 J2_system=systemmatrixunperturb*measurment_estimation_J2+J2_term*j2_acceleration;
 %calculate state from J2 model perspective
 J2_systemvalue=systemmatrixunperturb*J2_systemvalue+j2_acceleration;
 %calculate position and velocity error magnitude from J2 model
 J2_position(k-1)=norm(J2_systemvalue(1:3)-relative_J2LVLH);
 J2_velocity(k-1)=norm(J2_systemvalue(4:6)-gradient(relative_J2LVLH));
 %generate J2 model state with random noise 
 j2_model=[relative_J2LVLH-J2_position(k-1);gradient(relative_J2LVLH)-J2_velocity(k-1)];
 LOCAL_J2PVPH=(systemmatrixunperturb)*LOCAL_J2PVPH+(orbital_period/ceil(orbital_period/satrunspeed))*j2_acceleration+process_noiselinear;
 %position vector magnitude under j2
 J2_magnitude(k-1)=norm(LOCAL_J2PVPH(1:3));
 %velocity vector magnitude under j2
 J2_magnitudev(k-1)=norm(LOCAL_J2PVPH(4:6));
 %calculate jacobian
 deriviative_newstate=J2_system;
deriviative_newstate1=gradient(deriviative_newstate,measurment_estimation_J2(1));
deriviative_newstate2=gradient(deriviative_newstate,measurment_estimation_J2(2));
deriviative_newstate3=gradient(deriviative_newstate,measurment_estimation_J2(3));
deriviative_newstate4=gradient(deriviative_newstate,measurment_estimation_J2(4));
deriviative_newstate5=gradient(deriviative_newstate,measurment_estimation_J2(5));
deriviative_newstate6=gradient(deriviative_newstate,measurment_estimation_J2(6));
% % calculate whole system matrix
J2systemmatrix=transpose([deriviative_newstate1';deriviative_newstate2';deriviative_newstate3';deriviative_newstate4';deriviative_newstate5';deriviative_newstate6']);
%process noise for linear model
variacne_measurementJ2=1e8*[2.714e-6 3.7147e-6 4.714e-6 ];
R_matrix_linearJ2=variacne_measurementJ2(noise_sample)*eye(6);
%process noise covariance matrix
measurement_noiselinearJ2=1e10*diag([2.7e-7;2.7e-7;2.7e-7;2.7e-7;2.7e-7;2.7e-7]);
%measurement noise 
measurement_noise_linearJ2=sqrt(variacne_measurementJ2(noise_sample))*rand(1,6)';
%calculate measurement value
 y_linearJ2=linear_obervationmodel* LOCAL_J2PVPH+measurement_noise_linearJ2;
%initial estimation error and trace values of Predict error covariance
if k==2
error3=norm(measurment_estimation_J2(1:3)-LOCAL_J2PVPH(1:3));
total_errorj2(noise_sample,1)=real(20*log(error3/norm(LOCAL_J2PVPH(1:3))));
sumdiagonalj2(noise_sample,1)=trace(predicterror_cov_linearunperturb(1:3,1:3));
error33=norm(measurment_estimation_J2(4:6)-LOCAL_J2PVPH(4:6));
total_errorj21(noise_sample,1)=real(20*log(error33/norm(LOCAL_J2PVPH(4:6))));
sumdiagonalj21(noise_sample,1)=trace(predicterror_cov_linearunperturb(4:6,4:6));
 end 
 %predicted error covariance matrix estimate with J2
 predict_error_covJ2=J2systemmatrix*predicterror_cov_linearJ2*J2systemmatrix'+eye(6)*measurement_noiselinearJ2;
 estimate_measurmentJ2=J2systemmatrix*measurment_estimation_J2;
%calculate residual based on linear observation model
 residual_linearJ2=y_linearJ2-linear_obervationmodel*(estimate_measurmentJ2);
B=(linear_obervationmodel*predict_error_covJ2*transpose(linear_obervationmodel)+R_matrix_linearJ2);
 %calcualte kalman gain with J2 
kalman_gain_J2= predict_error_covJ2*transpose(linear_obervationmodel)*pinv(B);
%calculate residual norm of distacne and velocity portions
term=(eye(6)-kalman_gain_J2)*residual_linearJ2;
j2_residualdistance(k-1)=norm(term(1:3)+measurement_noise_linearJ2(1:3));
j2_residualvelocity(k-1)=norm(term(4:6)+measurement_noise_linearJ2(4:6));
j2_residual(k-1)=norm(residual_linearJ2(1:3));
j2_residual1(k-1)=norm(residual_linearJ2(4:6));
%difference between a prior state estimate and true state
a_postestimatej2=measurment_estimation_J2;
%calculate state estimation with J2 perturbation for linear model 
measurment_estimation_J2=estimate_measurmentJ2+kalman_gain_J2*residual_linearJ2;
j2_residualvar=LOCAL_J2PVPH-J2systemmatrix*a_postestimatej2;
%a posterior estimation error
j2var=(j2_model-a_postestimatej2);
%calculate predict error covariance matrix update for linear model
predicterror_cov_linearJ2=(eye(6)-kalman_gain_J2*linear_obervationmodel)*predict_error_covJ2;
 %TRACE of predict error covariance matrix
 if k>2
sumdiagonalj2(noise_sample,k-1)=trace(predicterror_cov_linearJ2(1:3,1:3));
%trace for velocity term
sumdiagonalj21(noise_sample,k-1)=trace(predicterror_cov_linearJ2(4:6,4:6));
 end
%% LVLH coordinates with J2
radial_track(k)=local_relativePOSITIONLVLH(1);
along_track(k)=local_relativePOSITIONLVLH(2);
cross_track(k)=local_relativePOSITIONLVLH(3);
%% update estimatin error array
%calculate estimation error with drag modified HCW model
if k>2
error=norm(measurment_estimation_linear(1:3)-relative_state_drag(1:3));
total_errordrag(noise_sample,k-1)=real(20*log(error/norm(FINAL_DRAG(1:3))));
%with respect to velocity
errordrag=norm(measurment_estimation_linear(4:6)-relative_state_drag(4:6));
total_errordrag1(noise_sample,k-1)=real(20*log(errordrag/norm(FINAL_DRAG(4:6))));
%calculate estimation error with normal HCW model relative_state
%with respect to poistion
error1=measurment_estimation_linearunperturb(1:3)- model_state(1:3);
total_errordragunperturb(noise_sample,k-1)=real(20*log(norm(error1)/norm(model_state(1:3))));
%with respect to velocity
error11=measurment_estimation_linearunperturb(4:6)- model_state(4:6);
total_errordragunperturb11(noise_sample,k-1)=real(20*log(norm(error11)/norm(model_state(4:6))));                                       
%calculate estimation error with J2 HCW modelrelative_state+final_coordinateECI_perturbationJ2LVLH
error3=norm(measurment_estimation_J2(1:3)-LOCAL_J2PVPH(1:3));
total_errorj2(noise_sample,k-1)=real(20*log(error3/norm(LOCAL_J2PVPH(1:3))));
%with respect to velocity
error33=norm(measurment_estimation_J2(4:6)-LOCAL_J2PVPH(4:6));
total_errorj21(noise_sample,k-1)=real(20*log(error33/norm(LOCAL_J2PVPH(4:6))));
% error2=norm(measurment_estimation_linearunperturb1(1:3)- LOCAL_J2PVPH(1:3));
% total_errordragunperturb1(noise_sample,k-1)=real(20*log(error3/norm(LOCAL_J2PVPH(1:3))))-10;
end

%%
    %update mean anomaly
     Mean_anomaly=Mean_anomaly+sqrt((GM/(a*a*a)))*satrunspeed;
 %% visualize deputy satellite motion for three cases
 %Deputy ECi coordinates
ECI_deputyunperturb=ECI_RELATIVE+([Xeci Yeci Zeci]');
ECI_deputyun1(k-1)=ECI_deputyunperturb(1);
ECI_deputyun2(k-1)=ECI_deputyunperturb(2);
ECI_deputyun3(k-1)=ECI_deputyunperturb(3);
ECI_deputydrag=-inv(DCM)*FINAL_DRAG(1:3)+([Xeci Yeci Zeci]');
ECI_deputyun11(k-1)=ECI_deputydrag(1);
ECI_deputyun22(k-1)=ECI_deputydrag(2);
ECI_deputyun33(k-1)=ECI_deputydrag(3);
ECI_deputyj2=inv(DCM)*LOCAL_J2PVPH(1:3)+([Xeci Yeci Zeci]');
ECI_deputyun111(k-1)=ECI_deputyj2(1);
ECI_deputyun222(k-1)=ECI_deputyj2(2);
ECI_deputyun333(k-1)=ECI_deputyj2(3);
 %Deputy LVLH coordinates
LVLH_deputyunperturb=DCM*ECI_RELATIVE;
LVLH_deputyun1(k-1)=LVLH_deputyunperturb(1);
LVLH_deputyun2(k-1)=LVLH_deputyunperturb(2);
LVLH_deputyun3(k-1)=LVLH_deputyunperturb(3);
%Deputy LVLH coordinates under drag model
LVLH_deputydrag=FINAL_DRAG(1:3);
LVLH_deputyun11(k-1)=LVLH_deputydrag(1);
LVLH_deputyun22(k-1)=LVLH_deputydrag(2);
LVLH_deputyun33(k-1)=LVLH_deputydrag(3);
%Deputy LVLH coordinates under j2 model
LVLH_deputyj2=LOCAL_J2PVPH(1:3);
LVLH_deputyun111(k-1)=LVLH_deputyj2(1);
LVLH_deputyun222(k-1)=LVLH_deputyj2(2);
LVLH_deputyun333(k-1)=LVLH_deputyj2(3);
 end
end
%% plot deputy satellite motion for three cases: unperturb drag force and J2 perturbation
figure(33)
subplot(2,1,1)
plot3 (ECI_deputyun1,ECI_deputyun2,ECI_deputyun3,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
hold on
plot3 (ECI_deputyun11,ECI_deputyun22,ECI_deputyun33,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','g','MarkerSize', 6);     
hold on
plot3 (ECI_deputyun111,ECI_deputyun222,ECI_deputyun333,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);     
xlabel('x-cor')
ylabel('y-cor')
zlabel('z-cor')
title('Unperturbed and Perturbed Satellite Motion in ECI frame')
legend('Deputy Satellite Motion with Drag Force','Deputy Unperturb Satellite Motion','Deputy Satellite Motion with J2 Perturbation')   
subplot(2,1,2)
plot3 (LVLH_deputyun1,LVLH_deputyun2,LVLH_deputyun3,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
hold on
plot3 (LVLH_deputyun11,LVLH_deputyun22,LVLH_deputyun33,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','g','MarkerSize', 6);     
hold on
plot3 (LVLH_deputyun111,LVLH_deputyun222,LVLH_deputyun333,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);     
xlabel('x-cor')
ylabel('y-cor')
zlabel('z-cor')
title('Unperturbed and Perturbed Satellite Motion in LVLH frame')
legend('Deputy Unperturb Satellite Motion','Deputy Satellite Motion with Drag Force','Deputy Satellite Motion with J2 Perturbation')   
%%
%average trace ratio for unperturbed model
b2=size(sumdiagonal);
b3=size(sumdiagonal1); 
NormRows1 = sum(sum(sumdiagonal,2))/(b2(1)*b2(2));
sumdiagonal = bsxfun(@rdivide,sumdiagonal,NormRows1);
NormRows = sum(sum(sumdiagonal1,2))/(b3(1)*b3(2));
sumdiagonal1 = bsxfun(@rdivide,sumdiagonal1,NormRows);
%%
%average trace ratio for drag-modified model
bdrag2=size(dragtrace);
bdrag3=size(dragtrace1);
NormRows11 = sum(sum(dragtrace,2))/(bdrag2(1)*bdrag2(2));
dragtrace = bsxfun(@rdivide,dragtrace,NormRows11);
NormRows22 = sum(sum(dragtrace1,2))/(bdrag3(1)*bdrag2(2));
dragtrace1 = bsxfun(@rdivide,dragtrace1,NormRows22);
%%
%normlize trace ratio for j2-modified model
j22=size(dragtrace);
j23=size(dragtrace1);
Norm = sum(sum(sumdiagonalj2,2))/(j22(1)*j22(2));
sumdiagonalj2 = bsxfun(@rdivide,sumdiagonalj2,abs(Norm));
NormRows222 = sum(sum(sumdiagonalj21,2))/(j23(1)*j23(2));
sumdiagonalj21 = bsxfun(@rdivide,sumdiagonalj21,abs(NormRows222));
%sum up position and velocity estimation error for every iteration 
final_total_errordragunperturb=(final_total_errordragunperturb+total_errordragunperturb);
final_total_errordragunperturb1=final_total_errordragunperturb1+total_errordragunperturb11;
final_total_errordrag=(final_total_errordrag+total_errordrag);
final_total_errordrag1=(final_total_errordrag1+total_errordrag1);
final_total_errorj2=(final_total_errorj2+total_errorj2);
final_total_errorj21=(final_total_errorj21+total_errorj21);
finalredisualdrag_norm=finalredisualdrag_norm+redisualdrag_norm(1,:);
finalredisualdrag_normvelocity=finalredisualdrag_normvelocity+redisualdrag_normvelocity(1,:);
finalredisualunper_normr=finalredisualunper_normr+redisual_normposition(1,:);
finalredisualunper_normv=finalredisualunper_normv+redisual_normvelocity(1,:);
finalredisualj2_normr=finalredisualj2_normr+j2_residual(1,:);
finalredisualj2_normv=finalredisualj2_normv+j2_residual1(1,:);
%sum up true position and velocity state norm for every iteration 
finallvlh_magnitude=lvlh_magnitude+finallvlh_magnitude;
finallvlh_magnitudev=lvlh_magnitudev+finallvlh_magnitudev;
finalunperturbpos=unperturbpos+finalunperturbpos;
finalunperturbvel=unperturbvel+finalunperturbvel;
end
%%averaged position and velocity estimation error and residual and true state norm for unperturbed
%%and perturbed case
total_errordragunperturb=final_total_errordragunperturb/local_simu_number;
total_errordragunperturb11=final_total_errordragunperturb1/local_simu_number;
total_errordrag=final_total_errordrag/local_simu_number;
total_errordrag1=final_total_errordrag1/local_simu_number;
total_errorj2=final_total_errorj2/local_simu_number;
total_errorj21=final_total_errorj21/local_simu_number;
redisualdrag_norm=finalredisualdrag_norm/local_simu_number;
redisualdrag_normvelocity=finalredisualdrag_normvelocity/local_simu_number;
redisual_normposition=finalredisualunper_normr/local_simu_number;
redisual_normvelocity=finalredisualunper_normv/local_simu_number;
j2_residual=finalredisualj2_normr/local_simu_number;
j2_residual1=finalredisualj2_normv/local_simu_number;
lvlh_magnitude=finallvlh_magnitude/local_simu_number;
lvlh_magnitudev=finallvlh_magnitudev/local_simu_number;
unperturbpos=finalunperturbpos/local_simu_number;
unperturbvel=finalunperturbvel/local_simu_number;


%% plot error pattern analysis under j2
figure(13)    
%subplot(2,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,j2_residual,'LineWidth',2);
hold on
plot(x,j2_residual1,'LineWidth',2);
hold off 
grid on
legend('||d_{r-k}||, unit:m','||d_{v-k}||,unit:m/s')
title('The Residual Norm for Position and Velocity:||d_{r-k}|| and ||d_{v-k}|| Comparison for J2-modified HCW Model')

%% plot drag error pattern analysis
figure(20) 
subplot(2,1,1)
x=linspace(0,length(redisual_normposition)-1,length(redisual_normposition));
x=x*orbital_period/length(x);
plot (x,drag_distance(1,:),'-.r*','LineWidth',2);
hold on
plot (x,drag_velocity(1,:),'-.r','LineWidth',2);
title('d_{r-metric} and d_{v-metric}')
legend('||d_{r-metric}||(unit:m)','||d_{v-metric}||(unit:m/s)')
xlabel('Time(unit:s)');
subplot(2,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,redisualdrag_norm(1,:),'LineWidth',2);
hold on
plot(x,redisualdrag_normvelocity(1,:),'LineWidth',2);
hold off 
grid on
legend('||d_{r-k}||,unit:m','||d_{v-k}||,unit:m/s')
title('The Residual Norm for Position and Velocity:||d_{r-k}|| and ||d_{v-k}|| Comparison for drag-modified HCW Model')
 %%   plot state norm comparison between drag model and unpuerturbed model
figure(2)    
subplot(2,1,1)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,lvlh_magnitudev,'LineWidth',2);
hold on
plot(x,drag_magnitudev,'LineWidth',2);
grid on
xlabel('Time(unit:m)'); 
legend('Unperturbed ||v_k||','Drag-perturbed ||v_k||')
ylabel('||v_k||(m/s)')
title('True Velocity Magnitude Comparison for Unperturbed and J2-perturbed Case')

xlabel('Time(unit:s)'); 
legend('Unperturbed ||r_k||','J2-perturbed ||r_k||')
subplot(2,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);

plot(x,lvlh_magnitude,'LineWidth',2);
hold on

plot(x,drag_magnitude,'LineWidth',2);
grid on
ylabel('||v_k||(m/s)')
title('True Position Magnitude Comparison for Unperturbed and Drag-perturbed Case')
xlabel('Time(unit:m)'); 
legend('Unperturbed ||r_k||','Drag-perturbed ||r_k||')
%% plot True state norm
figure(16)
subplot(2,1,1)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,lvlh_magnitude,'LineWidth',2);
grid on
title('True Position State Magnitude ||r_{k}||')
ylabel('||r_{k}||(m)')
xlabel('Time(s)')
legend('True Position State Norm')
subplot(2,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,lvlh_magnitudev,'LineWidth',2);
hold on
% plot(x,gradient(lvlh_magnitude-distance_error(1,:)),'LineWidth',2);
% hold off
grid on
xlabel('Time(s)')
title('True Velocity State Magnitude ||v_{k}|| ')
legend('True Velocity State Norm')
ylabel('||v_{k}||(m/s)')
%% plot state norm comparison between j2 model and unpuerturbed model
figure(3)
subplot(3,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,lvlh_magnitude,'LineWidth',2);
hold on
plot(x,J2_magnitude,'LineWidth',2);
grid on
ylabel('||v_k||(m/s)')
title('True Position Magnitude Comparison for Unperturbed and J2-perturbed Case')
xlabel('Time(unit:s)'); 
legend('Unperturbed ||r_k||','J2-perturbed ||r_k||')
subplot(3,1,1)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,lvlh_magnitudev,'LineWidth',2);
hold on
plot(x,J2_magnitudev,'LineWidth',2);
grid on
xlabel('Time(unit:s)'); 
legend('Unperturbed ||v_k||','J2-perturbed ||v_k||')
ylabel('||v_k||(m/s)')
title('True Velocity Magnitude Comparison for Unperturbed and J2-perturbed Case')
subplot(3,1,3)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,orbital_heightchief,'LineWidth',2);
hold on
plot(x,orbital_heightdeputy,'LineWidth',2);
grid on
xlabel('Time(unit:s)'); 
legend('||r_o|| for Chief Satellite','||r_o|| for Deputy Satellite')
title('Orbital Height Comparison for Chief and Deputy Satellite')
ylabel('||r_o||(m)')
%% plot ||d_{r-metric}|| and ||d_{v-metric}||
figure(5) 
subplot(2,1,1)
x=linspace(0,length(unperturbpos(1,:))-1,length(unperturbpos(1,:)));
x=x*orbital_period/length(x);
 plot (x,unperturbpos(1,:),'-.r*','LineWidth',2);
 hold on
plot (x,unperturbvel(1,:),'-.r','LineWidth',2);
title('||d_{r-metric}|| and ||d_{v-metric}||')
legend('||d_{r-metric}||(unit:m)','||d_{v-metric}||(unit:m/s)')
xlabel('Time(unit:s)');
subplot(2,1,2)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot(x,redisual_normposition,'LineWidth',2);
hold on
plot(x,redisual_normvelocity,'LineWidth',2);
hold off 
grid on
legend('||d_{r-k}||,unit:m','||d_{v-k}||,unit:m/s')
title('The Residual Norm for Position and Velocity:||d_{r-k}|| and ||d_{v-k}|| Comparison')
%% plot estimation error with unperturbed model
figure(9)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
subplot(3,1,1)
noise_color=['r','b','g','m','y','c','w','k','r'];
plot(x,total_errordragunperturb(1,:),'-.r*','LineWidth',2);
hold on
plot(x,total_errordragunperturb(2,:),'-.b*','LineWidth',2);
hold on
plot(x,total_errordragunperturb(3,:),'-.g*','LineWidth',2);
legend('\sigma_v^2=2.714e0 ','\sigma_v^2=2.714e1','\sigma_v^2=2.714e2')
ylabel('E_{r}');
title('Normalized unperturbed model Position Estimation Performance(\sigma_w^2=2.714e0)')
xlabel('Time(unit:s)');
subplot(3,1,2)
noise_color=['r','b','g','m','y','c','w','k','r']
plot(x,total_errordragunperturb11(1,:),'-.r','LineWidth',2);
hold on
plot(x,total_errordragunperturb11(2,:),'-.b','LineWidth',2);
hold on
plot(x,total_errordragunperturb11(3,:),'-.g','LineWidth',2);
legend('\sigma_v^2=2.714e0 ','\sigma_v^2=2.714e1','\sigma_v^2=2.714e2')
ylabel('E_{v}');
title('Normalized unperturbed model Velocity Estimation Performance(\sigma_w^2=2.714e0)')
xlabel('Time(unit:s)');
 subplot(3,1,3)
 x=linspace(0,k-2,k-1);
 x=x*orbital_period/length(x);
 plot(x,sumdiagonal(1,:),'-.r*','LineWidth',2)
hold on
plot(x,sumdiagonal(2,:),'-.b*','LineWidth',2)
hold on
plot(x,sumdiagonal(3,:),'-.g*','LineWidth',2)
hold on
plot(x,sumdiagonal1(1,:),'-.r','LineWidth',2)
hold on
plot(x,sumdiagonal1(2,:),'-.b','LineWidth',2)
hold on
plot(x,sumdiagonal1(3,:),'-.g','LineWidth',2)
title('trace(P_r) and trace(P_v)')
legend('trace(P_r)(\sigma_v^2=2.714e0)','trace(P_r)(\sigma_v^2=2.714e1)','trace(P_r)(\sigma_v^2=2.714e2)','trace(P_v)(\sigma_v^2=2.714e0)','trace(P_v)(\sigma_v^2=2.714e1)','trace(P_v)(\sigma_v^2=2.714e2)')
xlabel('Time(unit:s)');

 
%% plot estimation error with drag model
figure(10)
subplot(3,1,1)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
noise_color=['r','b','g','m','y','c','w','k','r'];
plot(x,total_errordrag(1,:),'-.r*','LineWidth',2);
hold on
plot(x,total_errordrag(2,:),'-.b*','LineWidth',2);
hold on
plot(x,total_errordrag(3,:),'-.g*','LineWidth',2);
grid on
% legend('\sigma_v^2=2.714e4','\sigma_v^2=3.714e4','\sigma_v^2=4.714e4')
legend('\sigma_v^2=2.714e4 ','\sigma_v^2=2.714e5','\sigma_v^2=2.714e6')
ylabel('E_{r}');
title('Normalized drag-modified model Deputy Satellite Position Estimation Performance(\sigma_w^2=2.714e2)')
xlabel('Time(unit:s)');
subplot(3,1,2)
 x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
noise_color=['r','b','g','m','y','c','w','k','r'];
plot(x,total_errordrag1(1,:),'-.r','LineWidth',2);
hold on
plot(x,total_errordrag1(2,:),'-.b','LineWidth',2);
hold on
plot(x,total_errordrag1(3,:),'-.g','LineWidth',2);
grid on
legend('\sigma_v^2=2.714e4','\sigma_v^2=2.714e5','\sigma_v^2=2.714e6')
ylabel('E_{v}');
title('Normalized drag-modified model Deputy Satellite Velocity Estimation Performance(\sigma_w^2=2.714e2)')
xlabel('Time(unit:s)');
subplot(3,1,3)
x=linspace(0,k-2,k-1);
 x=x*orbital_period/length(x);
 plot (x,dragtrace(1,:),'-.r*','LineWidth',2);
 hold on
plot (x,dragtrace(2,:),'-.b*','LineWidth',2);
hold on
plot (x,dragtrace(3,:),'-.g*','LineWidth',2);
hold on
plot (x,dragtrace1(1,:),'-.r','LineWidth',2);
hold on
plot (x,dragtrace1(2,:),'-.b','LineWidth',2);
hold on
plot (x,dragtrace1(3,:),'-.g','LineWidth',2);
%  legend('\sigma_w^2(Relative Distance Error)) ','\sigma_w^2(Relative Distance Error) ','\sigma_w^2(Relative Distance Error) ')  
% legend('Position(\sigma_v^2=2.714e1)','Position(\sigma_v^2=3.714e1)','Position(\sigma_v^2=4.714e1)','Velocity(\sigma_v^2=2.714e1)','Velocity(\sigma_v^2=3.714e1)','Velocity(\sigma_v^2=4.714e1)')
legend('trace(P_r)(\sigma_v^2=2.714e4)','trace(P_r)(\sigma_v^2=2.714e5)','trace(P_r)(\sigma_v^2=2.714e6)','trace(P_v)(\sigma_v^2=2.714e4)','trace(P_v)(\sigma_v^2=2.714e5)','trace(P_v)(\sigma_v^2=2.714e6)')
title('trace(P_r) and trace(P_v)')
xlabel('Time(unit:s)');
%% plot j2 modified model performance
figure(11)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
noise_color=['r','b','g','m','y','c','w','k','r']
subplot(3,1,1)
plot(x(1:k-1),total_errorj2(1,1:k-1),'-.ro','LineWidth',2);
hold on
plot(x(1:k-1),total_errorj2(2,1:k-1),'-.bo','LineWidth',2);
hold on
plot(x(1:k-1),total_errorj2(3,1:k-1),'-.go','LineWidth',2);
grid on
hold off
legend('\sigma_v^2=2.714e2','\sigma_v^2=3.714e2','\sigma_v^2=4.714e2')
ylabel('E_{r}');
title('Normalized Position Estimation Performance with J2 modified Model(\sigma_w^2=2.714e3')
subplot(3,1,2)
plot(x(1:k-1),total_errorj21(1,1:k-1),'-.ro','LineWidth',2);
hold on
plot(x(1:k-1),total_errorj21(2,1:k-1),'-.bo','LineWidth',2);
hold on
plot(x(1:k-1),total_errorj21(3,1:k-1),'-.go','LineWidth',2);
grid on
hold off
legend('\sigma_v^2=2.714e2','\sigma_v^2=3.714e2','\sigma_v^2=4.714e2')
ylabel('E_{v}');
title('Normalized Velocity Estimation Performance with J2 modified Model(\sigma_w^2=2.714e3)')
xlabel('Time(unit:s)');
subplot(3,1,3)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);
plot (x,sumdiagonalj21(1,:),'-.r*','LineWidth',2);
hold on
plot (x,sumdiagonalj21(2,:),'-.b*','LineWidth',2);
hold on
plot (x,sumdiagonalj21(3,:),'-.g*','LineWidth',2);
hold on
plot (x,sumdiagonalj2(1,:),'-.r','LineWidth',2);
hold on
plot (x,sumdiagonalj2(2,:),'-.b','LineWidth',2);
hold on
plot (x,sumdiagonalj2(3,:),'-.g','LineWidth',2);
grid on
legend('trace(P_r)(\sigma_v^2=2.714e2)','trace(P_r)(\sigma_v^2=3.714e2)','trace(P_r)(\sigma_v^2=4.714e2)','trace(P_v)(\sigma_v^2=2.714e2)','trace(P_v)(\sigma_v^2=3.714e2)','trace(P_v)(\sigma_v^2=4.714e2)')
title('trace(P_r) and trace(P_v)')
xlabel('Time(unit:s)');
end