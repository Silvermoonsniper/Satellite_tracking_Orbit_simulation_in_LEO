%% This function aims to generate a circular orbit with specific initialization and return corresponding orbital elements
%input parameters:
%initial velocity unit vector and VELOCITY magnitude
function [a,inc,RAAN,w,M0,e,R_magnitude,v,r]=orbital_elements(v,vmag)
% Planetary gravitational constant for Earth, (mu = GMearth) (m^3/s^2)
mu = 398.6004418e12;  
R_magnitude=mu/vmag^2; 
r=[R_magnitude 0 0]; 
v= vmag*(v/sqrt(dot(v, v))); %Velocity Vector for sat 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERTING STATE VECTORS INTO ORBITAL ELEMENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%=
V=v;
rmag = sqrt(dot(r, r)); %Position Magnitude
vmag = sqrt(dot(v, v)); %Velocity Magnitude

rhat = r/rmag; %Position Unit Vector
vhat = v/vmag; %Velocity Unit Vector

hv = cross(r, v); %Angular Momentum Vector
hmag = sqrt(dot(hv, hv)); %Angular Momentum Magnitude
hhat = hv/hmag; %Angular Momentum Unit Vector

%Eccentricity Vector
vtmp = v / mu;
ecc = cross(vtmp, hv);
ecc = ecc - rhat;
p=hmag^2/mu;
%SEMIMAJOR AXIS (a)
a = 1 / (2 / rmag - vmag * vmag / mu);
%ECCENTRICITY (e) %0<e<1
e = sqrt(1-p/a);
%calculate eccentric anomaly
E=acos((1-R_magnitude/a)/e);
%INCLINATION (inc) %in rad
inc = acos(hhat(3));


%RIGHT ASCENSION OF ASCENDING NODE (RAAN) %in rad
if (inc > 0.00000001)
    RAAN = asin(hhat(1)/sin(inc));
else
   RAAN = 0;
end
%ARGUMENT OF PERIGEE (w) %in rad
if (e > 0.00000001)
   w=cross([0  0 1],h); 
else
   w = 0;
end
%MEAN ANOMALY (M0)
M0 =E-e*sin(E); %in rad
end