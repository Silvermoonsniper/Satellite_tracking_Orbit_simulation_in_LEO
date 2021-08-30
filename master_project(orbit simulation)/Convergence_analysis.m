%% This function is aiming to invertigate how process noise influence the predict error variation in the 
%% extended kalman filter we apply for satellite tracking with presence of J2 perturbation with a fixed 
%% range of observation noise variances
function Convergence_analysis
%% initialization of deputy satellite in LVLH frame
v1=[0 0 1];
vmag=7880;
%revolution number
rev_number=1;
%satellite runspeed
satrunspeed=300;
%earth radius
EARTH_RAIDUS=6371000;% unit:m
%Gravitational Constant
GM=398.6004418e12;
[a,Inc,omega,w,Mean_anomaly,e,R_magnitude,v1,r1]=orbital_elements(v1,vmag);
%angular velocity of chief satellite
angular_v=sqrt(GM/norm(R_magnitude)^3);
%initial position vector in LVLH frame
delta_r=[1e3/6 0 1e3];
%initial velocity vector in LVLH frame
delta_v=[0 -2*angular_v*delta_r(1) 0];
%calculate parameter for HCW linearized solution
phase_x=atan(real((1/angular_v)*(-delta_v(1))/(3*delta_r(1))));
phase_z=atan(-(1/angular_v)*(delta_v(3)/delta_r(3)));
A_z=delta_r(3)/real(cos((phase_z)));
A_x=-(3*delta_r(1)+2*delta_r(2)/angular_v)/cos(real(phase_x));
orbital_period=sqrt((4*pi*pi*R_magnitude^3)/GM);
variacne_processJ2=[2.714e-1 2.714e0 2.714e1 2.714e2 2.714e3];
for noise_sample1=1:5
  
k=1;
%initial state guess and predict error covariance matrix for J2 model
ww=1e-6*[1.8e6, 1.8e6, 1.8e6, 18 ,18 , 18];
measurment_estimation_J2=rand([1 6])';
predicterror_cov_linearJ2=diag(ww);
LOCAL_J2PVPH=[delta_r';delta_v'];
for time_instant=1:rev_number*ceil(orbital_period/satrunspeed)
    
    k=k+1;
 %calculate angular velocity of both satellites same for both due to same
 %semi-major axis
  angular_v=sqrt(GM/R_magnitude^3);
  %calculate eccentric anomaly from mean anomaly
    E_new=real(Mean_anomaly);
    for i=1:5
         E_new=Mean_anomaly+(Mean_anomaly + e*sin(E_new) - E_new)/(1 - e*cos(E_new));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    %true anomaly
     true_anomaly=2*atan(sqrt((1+e)/(1-e))*tan(E_new/2));
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

%random noise to generate true state
process_noiselinear=sqrt(2.7e1)*rand(1,6)';
%establish system matrix f(orbital_period/ceil(orbital_period/satrunspeed))*
 systemmatrixunperturb=([zero_matrix Identity;
     k_matrix+J  omeag_matrix ])+eye(6);

%% J2 perturbed extended kalman filter
 %J2 acceleration components in LVLH frame
 J2=1.083e-3;
 j2_factor=3*GM*J2*EARTH_RAIDUS^2/(2*R_magnitude^4);
 %angle theta
 theta=w+true_anomaly;
 %additional acceleration caused by J2
 j2_acceleration=j2_factor*[0;0;0;-1+3*sin(Inc)^2*sin(theta)^2;-2*sin(Inc)^2*sin(theta)*cos(theta);-2*sin(Inc)*cos(theta)*sin(theta)];
 
 %J2 distance variation along three directions solved by J2 modified HCW
 %equation
J2_term=3*J2*EARTH_RAIDUS^2/(2*R_magnitude);
J2_radial(k-1)=J2_term*((1/3)*sin(Inc)^2*cos(theta)^2+((1/3)*sin(Inc)^2-1)+(1-(2/3)*sin(Inc)^2)*cos(theta));
J2_ALONG(k-1)=J2_term*((2/3)*sin(Inc)^2*cos(theta)+(1-(2/3)*sin(Inc)^2)*cos(theta))+(j2_factor*sin(Inc)^2*cos(2*theta))/(4*angular_v^2);
J2_cross(k-1)=j2_factor*sin(Inc)*sin(2*theta)/(-5*angular_v^2);
%true position vector under J2(from physical model) 
relative_J2LVLH=[J2_radial(k-1);J2_ALONG(k-1);J2_cross(k-1)];
 %modified J2 perturbed HCW model

 J2_system=systemmatrixunperturb*measurment_estimation_J2+J2_term*j2_acceleration;
 %calculate state from J2 model perspective
%  J2_systemvalue=systemmatrixunperturb*J2_systemvalue+j2_acceleration;
 %generate J2 model state with random noise 
  LOCAL_J2PVPH=(systemmatrixunperturb)*LOCAL_J2PVPH+(orbital_period/ceil(orbital_period/satrunspeed))*j2_acceleration+process_noiselinear;
deriviative_newstate=J2_system;
deriviative_newstate1=gradient(deriviative_newstate,measurment_estimation_J2(1));
deriviative_newstate2=gradient(deriviative_newstate,measurment_estimation_J2(2));
deriviative_newstate3=gradient(deriviative_newstate,measurment_estimation_J2(3));
deriviative_newstate4=gradient(deriviative_newstate,measurment_estimation_J2(4));
deriviative_newstate5=gradient(deriviative_newstate,measurment_estimation_J2(5));
deriviative_newstate6=gradient(deriviative_newstate,measurment_estimation_J2(6));
%linear_obervationmodel
linear_obervationmodel=eye(6);
% % calculate whole system matrix
J2systemmatrix=transpose([deriviative_newstate1';deriviative_newstate2';deriviative_newstate3';deriviative_newstate4';deriviative_newstate5';deriviative_newstate6']);
%process noise for linear model
variacne_measurementJ2=1e8*[2.714e-6 3.7147e-6 4.714e-6 ];
noise_sample=1;
R_matrix_linearJ2=variacne_measurementJ2(noise_sample)*eye(6);
%process noise covariance matrix
measurement_noiselinearJ2=variacne_processJ2(noise_sample1)*eye(6);
% measurement noise 
measurement_noise_linearJ2=sqrt(variacne_measurementJ2(noise_sample))*rand(1,6)';
%calculate measurement value
 y_linearJ2=linear_obervationmodel* LOCAL_J2PVPH+measurement_noise_linearJ2;


 sumdiagonalj2(noise_sample1,k-1)=trace(predicterror_cov_linearJ2(1:3,1:3)); 
 sumdiagonalj21(noise_sample1,k-1)=trace(predicterror_cov_linearJ2(4:6,4:6)); 
 %predicted error covariance matrix estimate for linear modelmeasurement_noiselinearJ2
 predict_error_covJ2=J2systemmatrix*predicterror_cov_linearJ2*J2systemmatrix'+eye(6)*measurement_noiselinearJ2;
 estimate_measurmentJ2=J2systemmatrix*measurment_estimation_J2;
%calculate residual based on linear observation model
 residual_linearJ2=y_linearJ2-linear_obervationmodel*(estimate_measurmentJ2);
B=(linear_obervationmodel*predict_error_covJ2*transpose(linear_obervationmodel)+R_matrix_linearJ2);
 %calcualte kalman gain
kalman_gain_J2= predict_error_covJ2*transpose(linear_obervationmodel)*pinv(B);
%calculate residual norm of distacne and velocity portions
term=(eye(6)-kalman_gain_J2)*residual_linearJ2;

%calculate state estimation with J2 perturbation for linear model 
measurment_estimation_J2=estimate_measurmentJ2+kalman_gain_J2*residual_linearJ2;
%calculate predict error covariance matrix update for linear model
predicterror_cov_linearJ2=(eye(6)-kalman_gain_J2*linear_obervationmodel)*predict_error_covJ2;
%update mean anomaly
Mean_anomaly=Mean_anomaly+sqrt((GM/(a*a*a)))*satrunspeed;


end
end
%plot trace results
figure(1)
x=linspace(0,k-2,k-1);
x=x*orbital_period/length(x);

 plot (x,sumdiagonalj21(1,:),'-.r*','LineWidth',2);
 hold on
plot (x,sumdiagonalj21(2,:),'-.b*','LineWidth',2);
hold on
plot (x,sumdiagonalj21(3,:),'-.g*','LineWidth',2);
hold on
plot (x,sumdiagonalj21(4,:),'-.c*','LineWidth',2);
 hold on
plot (x,sumdiagonalj21(5,:),'-.y*','LineWidth',2);
hold on
plot (x,sumdiagonalj2(1,:),'-.r','LineWidth',2);
 hold on
plot (x,sumdiagonalj2(2,:),'-.b','LineWidth',2);
hold on
plot (x,sumdiagonalj2(3,:),'-.g','LineWidth',2);
hold on
plot (x,sumdiagonalj2(4,:),'-.c','LineWidth',2);
 hold on
plot (x,sumdiagonalj2(5,:),'-.y','LineWidth',2);
grid on
legend('trace(P_r)(\sigma_w^2=2.714e-1)','trace(P_r)(\sigma_w^2=2.714e0)','trace(P_r)(\sigma_w^2=2.714e1)','trace(P_r)(\sigma_w^2=2.714e2)','trace(P_r)(\sigma_w^2=2.714e3)','trace(P_v)(\sigma_w^2=2.714e-1)','trace(P_v)(\sigma_w^2=2.714e0)','trace(P_v)(\sigma_w^2=2.714e1)','trace(P_v)(\sigma_w^2=2.714e2)','trace(P_v)(\sigma_w^2=2.714e3)')
 
xlabel('Time(unit:s)');
title('trace(P_r) and trace(P_v)')

end