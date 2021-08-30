%% This function aims to simulate satellite orbit in LEO and investigate how velocity magnitude influence 
%% satellite orbital height variation, also satellite motion in ECEF frame is also done here
function  orbit(x_intial,y_intial,z_intial)
%intial position vector
r=[x_intial y_intial z_intial];
r1=[x_intial y_intial z_intial];
%intial velocity vector
v=[0 1 1];
v1=[0 1 1];
v2=[0 1 1];
v3=[0 1 1];
%intial velocity magnitude
% vmag=5000 %unit:m/s
vmag=7800;
vmag1=7850;
vmag2=7900;
vmag3=7950;
%intial velocity vector
v=v*(vmag/sqrt(dot(v,v)));
v1=v1*(vmag1/sqrt(dot(v1,v1)));
v2=v2*(vmag2/sqrt(dot(v2,v2)));
v3=v3*(vmag3/sqrt(dot(v3,v3)));

%Number of Days since Jan 1, 2000
J2000_days=103752/24;
%satellite speed
satrunspeed=200; 
%revolution time
rev_number=1;
r=r*1000;
r1=r1*1000;
%Gravitational Constant
GM=398.6004418e12;
%earth radius
EARTH_RAIDUS=6371000;% unit:m
% the semi-major axis can be derived from keplerian 3rd Law
R_magnitude=sqrt(dot(r,r));
%semi_major axis:a
a=1/((2/R_magnitude)-vmag^2/GM);
a1=1/((2/R_magnitude)-vmag1^2/GM);
a2=1/((2/R_magnitude)-vmag2^2/GM);
a3=1/((2/R_magnitude)-vmag3^2/GM);
rhat = r/R_magnitude; %Position Unit Vector
vhat = v/vmag; %Velocity Unit Vector
% angular momentum vector and magnitude
h=cross(r,v);
h_magnitude=sqrt(dot(h,h));
h1=cross(r,v1);
h_magnitude1=sqrt(dot(h1,h1));
h2=cross(r,v2);
h_magnitude2=sqrt(dot(h2,h2));
h3=cross(r,v3);
h_magnitude3=sqrt(dot(h3,h3));
h_unitvector=h/h_magnitude;
h_unitvector1=h1/h_magnitude1;
h_unitvector2=h2/h_magnitude2;
h_unitvector3=h3/h_magnitude3;

% eccentricity vector
p=h_magnitude^2/GM;
ECC=cross(v,h)/GM;
ECC=ECC-r;
p1=h_magnitude1^2/GM;
ECC1=cross(v1,h1)/GM;
ECC1=ECC1-r;
p2=h_magnitude2^2/GM;
ECC2=cross(v2,h2)/GM;
ECC2=ECC2-r;
p3=h_magnitude3^2/GM;
ECC3=cross(v3,h3)/GM;
ECC3=ECC3-r;
%eccentricity 
e=sqrt(1-p/a);
e1=sqrt(1-p1/a1);
e2=sqrt(1-p2/a2);
e3=sqrt(1-p3/a3);
%calculate eccentric anomaly
E=acos((1-R_magnitude/a)/e);
E1=acos((1-R_magnitude/a)/e1);
E2=acos((1-R_magnitude/a)/e2);
E3=acos((1-R_magnitude/a)/e3);
%ECCENTRICITY
% e=abs(((p/R_magnitude)-1)/cos(true_anomaly));
%calculate local velocity of satellite due to angular momentum conservation
total_energy=0.5*vmag^2-GM/(R_magnitude*R_magnitude);

% mean anomaly
Mean_anomaly=E-e*sin(E);
Mean_anomaly1=E1-e*sin(E1);
Mean_anomaly2=E2-e*sin(E2);
Mean_anomaly3=E3-e*sin(E3);
% inclination angle
 Inc=acos(h(3)/h_magnitude);
 Inc1=acos(h1(3)/h_magnitude1);
 Inc2=acos(h2(3)/h_magnitude2);
 Inc3=acos(h3(3)/h_magnitude3);
% ascending node at longtitude
omega=asin(h_unitvector(1)/sin(Inc));
omega1=asin(h_unitvector1(1)/sin(Inc1));
omega2=asin(h_unitvector2(1)/sin(Inc2));
omega3=asin(h_unitvector3(1)/sin(Inc3));
%  argument of perigee W
n=cross([0  0 1],h);
n1=cross([0  0 1],h1);
n2=cross([0  0 1],h2);
n3=cross([0  0 1],h3);
%magnitude of n vector
n_mag=sqrt(dot(n,n));
n_mag1=sqrt(dot(n1,n1));
n_mag2=sqrt(dot(n2,n2));
n_mag3=sqrt(dot(n3,n3));

% evector=(1/GM)*((vmag^2-GM/R_magnitude)*r-dot(r,v)*v);
evector=(cross(v,h)/GM)-r/R_magnitude;
evector1=(cross(v1,h1)/GM)-r/R_magnitude;
evector2=(cross(v2,h2)/GM)-r/R_magnitude;
evector3=(cross(v3,h3)/GM)-r/R_magnitude;
%magnitude of eccentricity vector
emag=sqrt(dot(evector,evector));
emag1=sqrt(dot(evector1,evector1));
emag2=sqrt(dot(evector2,evector2));
emag3=sqrt(dot(evector3,evector3));
 w=angle(acos(dot(n,evector)/emag*n_mag));
 w1=angle(acos(dot(n1,evector1)/emag1*n_mag1));
 w2=angle(acos(dot(n2,evector2)/emag2*n_mag2));
 w3=angle(acos(dot(n3,evector3)/emag3*n_mag3));
%w=1.2;
%  w=pi-omega;
%adjust argument of perigee
if evector(3)<0
    w=2*pi-w;
end
if evector1(3)<0
    w1=2*pi-w1;
end
if evector2(3)<0
    w2=2*pi-w2;
end
if evector3(3)<0
    w3=2*pi-w3;
end
 
%ERROR if Launch position is inside the Earth
if (sqrt (r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) <= 6731000)
   % blast(r(1), r(2), r(3), 2000000);
    ErrorMsg='Launch Position Inside Earth'
    return;
end
%check eccentricity validity
if e>1 || e<0
    ErrorMsg='eccentricity invalid';
end
%modify mean anomaly
if (Mean_anomaly>2*pi)
    Mean_anomaly=Mean_anomaly-2*pi;
end
if (Mean_anomaly1>2*pi)
    Mean_anomaly1=Mean_anomaly1-2*pi;
end
if (Mean_anomaly2>2*pi)
    Mean_anomaly2=Mean_anomaly2-2*pi;
end
if (Mean_anomaly3>2*pi)
    Mean_anomaly3=Mean_anomaly3-2*pi;
end

%Plotting the ECI Axes
lim=(1+e)*a;
clf
axis([-lim, lim, -lim, lim, -lim, lim])	
view(150,15) 
axis equal
shg
hold on
grid on
title('Satellite Orbital Simulation with Different Velocity Magnitudes');
line([0 lim],[0 0],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 lim],[0 0],'Color', 'black', 'Marker','.','LineStyle','-')
line([0 0],[0 0],[0 lim],'Color', 'black', 'Marker','.','LineStyle','-')
%plot the intial velocity and position vector
line([r(1)+v(1)*1000],[r(2)+v(2)*1000],[r(3)+v(3)*1000]);
line(r(1),r(2),r(3));

% Plotting the Earth
equat_rad=6378137.00;
polar_rad=6356752.3142;
[xx yy zz]=ellipsoid (0,0,0,equat_rad, equat_rad, polar_rad);
load('topo.mat','topo','topomap1');
    topo2 = [topo(:,181:360) topo(:,1:180)];
    pro.FaceColor= 'texture';
    pro.EdgeColor = 'none';
    pro.FaceLighting = 'phong';
    pro.Cdata = topo2;
   earth= surface(xx,yy,zz,pro);
    colormap(topomap1)
omega_earth = 7.292115855377074e-005; % (rad/sec)  
Go = 1.727564365843028; % (rad)  
GMST = Go + omega_earth*86400*(J2000_days + 0.5);
GMST = GMST - 2*pi*floor(GMST/(2*pi));
GMST_deg=GMST*(180/pi);

% draw dynamic satellite motion
orbital_period=sqrt((4*pi*pi*a*a*a)/GM);
%represent intial position coordinate and radius value
x_cor(1)=r(1);
y_cor(2)=r(2);
z_cor(3)=r(3);
x_cor1(1)=r(1);
y_cor1(2)=r(2);
z_cor1(3)=r(3);
x_cor2(1)=r(1);
y_cor2(2)=r(2);
z_cor2(3)=r(3);
x_cor3(1)=r(1);
y_cor3(2)=r(2);
z_cor3(3)=r(3);
status(1)=1;

r(1)=R_magnitude;
k=1;
%intial relative position
 relative_d(1)=0;
   relative_d1(1)=0;
   relative_d2(1)=0;
   %relative coordinate representation
   relativex1(1)=0;
   relativex2(1)=0;
   relativex3(1)=0;
   relativey1(1)=0;
   relativey2(1)=0;
   relativey3(1)=0;
   relativez1(1)=0;
   relativez2(1)=0;
   relativez3(1)=0;
%plot the intial relative velocity and position vector
% line(relativex1(1)+v(1)*1000],[r(2)+v(2)*1000],[r(3)+v(3)*1000]);
line(relativex1(1),relativey1(1),relativez1(1));
   %calculate iteration number
iterationnumber=rev_number*ceil(orbital_period/satrunspeed);
%bulid an empty array to place velocity values
v_local(1)=vmag;
 %intial minimum speed of satellite at given orbit height
  minium_speed(1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r(1))));
  minium_speed1(1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r(1))));
  minium_speed2(1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r(1))));
  minium_speed3(1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r(1))));
for time_instant=1:rev_number*ceil(orbital_period/satrunspeed)
    
    k=k+1;
     
    %calculate eccentric anomaly from mean anomaly
    E_new=Mean_anomaly;
    for i=1:5
         E_new=Mean_anomaly+(Mean_anomaly + e*sin(E_new) - E_new)/(1 - e*cos(E_new));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    E_new1=Mean_anomaly1;
    for i=1:5
         E_new1=Mean_anomaly1+(Mean_anomaly1 + e*sin(E_new1) - E_new1)/(1 - e*cos(E_new1));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    E_new2=Mean_anomaly2;
    for i=1:5
         E_new2=Mean_anomaly2+(Mean_anomaly2 + e*sin(E_new2) - E_new2)/(1 - e*cos(E_new2));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    E_new3=Mean_anomaly3;
    for i=1:5
         E_new3=Mean_anomaly3+(Mean_anomaly3 + e*sin(E_new3) - E_new3)/(1 - e*cos(E_new3));
%         E_new=Mean_anomaly+e*sin(E_new);
    end
    %calculate true anomaly
    true_anomaly=2*atan(sqrt((1+e)/(1-e))*tan(E_new/2));
    %compute radius vector norm ||r|| 
    R=a*(1-e*sin(E_new));
    r_magnitude(k-1)=a*(1-e*sin(E_new));
    R1=a1*(1-e1*sin(E_new1));
    r1_magnitude(k-1)=a*(1-e1*sin(E_new));
    R2=a2*(1-e*sin(E_new));
    r2_magnitude(k-1)=a*(1-e2*sin(E_new));
    R3=a3*(1-e3*sin(E_new3));
    r3_magnitude(k-1)=a*(1-e3*sin(E_new));
   %check whether satellite is dropping down to earth
    if r_magnitude(k-1)<=EARTH_RAIDUS
       error='satellite drops'
       return;
   end
   if r1_magnitude(k-1)<=EARTH_RAIDUS
       error='satellite drops'
       return;
   end
   if r2_magnitude(k-1)<=EARTH_RAIDUS
       error='satellite drops'
       return;
   end
   if r3_magnitude(k-1)<=EARTH_RAIDUS
       error='satellite drops'
       return;
   end
   %construct rotation matrix to convert ECI coordinates into ECEF
    %coordinatesGMST 
    angularvelocity_inradians= 7.2921e-5;
    rotation_ECItoECEF=[cos(-GMST-angularvelocity_inradians*k) sin(-GMST-angularvelocity_inradians*k) 0;
     -sin(-GMST-angularvelocity_inradians*k) cos(-GMST-angularvelocity_inradians*k) 0;
     0 0 1];
   
   %satellite orbit height calculation
    satelliteheight(k-1)=r_magnitude(k-1)-EARTH_RAIDUS;
    satelliteheight1(k-1)=r1_magnitude(k-1)-EARTH_RAIDUS;
    satelliteheight2(k-1)=r2_magnitude(k-1)-EARTH_RAIDUS;
    satelliteheight3(k-1)=r3_magnitude(k-1)-EARTH_RAIDUS;
  %minimum speed of satellite at given orbit height
  minium_speed(k-1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r_magnitude(k-1))));
  minium_speed1(k-1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r1_magnitude(k-1))));
  minium_speed2(k-1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r2_magnitude(k-1))));
  minium_speed3(k-1)=sqrt(2*GM*((1/EARTH_RAIDUS)-(1/r3_magnitude(k-1))));
    Xeci = R*(cos(w + true_anomaly)*cos(omega) - sin(w+true_anomaly)*sin(omega)*cos(Inc));
    Yeci= R*(cos(w + true_anomaly)*sin(omega) + sin(w+true_anomaly)*cos(omega)*cos(Inc));
    Zeci = R*(sin(w + true_anomaly)*sin(Inc));
 %ECEF coordinates
    coordinate_ECEF= rotation_ECItoECEF*[Xeci;Yeci;Zeci];
    
    Xeci1 = R1*(cos(w1 + true_anomaly)*cos(omega1) - sin(w1+true_anomaly)*sin(omega1)*cos(Inc1));
    Yeci1= R1*(cos(w1 + true_anomaly)*sin(omega1) + sin(w1+true_anomaly)*cos(omega1)*cos(Inc1));
    Zeci1 = R1*(sin(w1 + true_anomaly)*sin(Inc1));
     %ECEF coordinates
          coordinate_ECEF1= rotation_ECItoECEF*[Xeci1;Xeci1;Xeci1];
     Xeci2 = R2*(cos(w2 + true_anomaly)*cos(omega2) - sin(w2+true_anomaly)*sin(omega2)*cos(Inc2));
    Yeci2= R2*(cos(w2 + true_anomaly)*sin(omega2) + sin(w2+true_anomaly)*cos(omega2)*cos(Inc2));
    Zeci2 = R2*(sin(w2 + true_anomaly)*sin(Inc2));
     %ECEF coordinates
          coordinate_ECEF2= rotation_ECItoECEF*[Xeci2;Yeci2;Zeci2];
     Xeci3 = R3*(cos(w3 + true_anomaly)*cos(omega3) - sin(w3+true_anomaly)*sin(omega3)*cos(Inc3));
    Yeci3= R3*(cos(w3 + true_anomaly)*sin(omega3) + sin(w3+true_anomaly)*cos(omega3)*cos(Inc3));
    Zeci3 = R3*(sin(w3 + true_anomaly)*sin(Inc3));
     %ECEF coordinates
          coordinate_ECEF3= rotation_ECItoECEF*[Xeci3;Xeci3;Xeci3];
    %raidusmagnitude=sqrt(Xeci*Xeci+Yeci*Yeci+Zeci*Zeci);
%    figure(1)
   %relative distance calculation of satellite 1 with respect to others
   relative_d(k)=sqrt((Xeci-Xeci1)*(Xeci-Xeci1)+(Yeci-Yeci1)*(Yeci-Yeci1)+(Zeci-Zeci1)*(Zeci-Zeci1));
   relative_d1(k)=sqrt((Xeci-Xeci2)*(Xeci-Xeci2)+(Yeci-Yeci2)*(Yeci-Yeci2)+(Zeci-Zeci2)*(Zeci-Zeci2));
   relative_d2(k)=sqrt((Xeci-Xeci3)*(Xeci-Xeci3)+(Yeci-Yeci3)*(Yeci-Yeci3)+(Zeci-Zeci3)*(Zeci-Zeci3));
   %relative coordinate representation
   relativex1(k)=Xeci1-Xeci;
   relativex2(k)=Xeci2-Xeci;
   relativex3(k)=Xeci3-Xeci;
   relativey1(k)=Yeci1-Yeci;
   relativey2(k)=Yeci2-Yeci;
   relativey3(k)=Yeci3-Yeci;
   relativez1(k)=Zeci1-Zeci;
   relativez2(k)=Zeci2-Zeci;
   relativez3(k)=Zeci3-Zeci;
%% this commented part investigate controlability in Lyapunove sense, not part of project
% orbital stability analysis
%   x1=r(k);
%   x2= sqrt(v_local(k)^2-GM*p/r(k)^2);
%   x4=sqrt(GM*p)/r(k)^2;
%   x3=Mean_anomaly;
%    input_matrix=[0  1  0  0;
%        x4^2+2*GM/x1^3  0  0  2*x1*x4;
%        0  0 0  1;
%        2*x2*x4/x1^2  -2*x4/x1  0  -2*x2/x1];
%    %calculate characteristic equation
%    
%    charcteristic_poly=charpoly( input_matrix);
%    %calculate eigenvalue
%    eigenvalue=eig(input_matrix);
%    
%    %calculate real part of eigenvalues
%    realpart=real(eigenvalue);
%    number=0;
   %if all eigenvalues contain positive real part, then it is unstable
%    for i=1:4
%        if realpart(i)>0
%            number=number+1;
%        else
%            number=number;
%        end
%    end
%    if number==1
%        
%        stabilty='orbit is unstable'
%        status(k)=0;
%    else
%        stabilty='orbit is stable'
%         status(k)=1;
%    end
  %pause for 0.01s to continue iteration
    pause(0.1);
% append ECI and ECEF coordinates
    if time_instant>=1 && time_instant<ceil(orbital_period/satrunspeed+1) 
          x_cor(k)=Xeci;
          y_cor(k)=Yeci;
          z_cor(k)=Zeci;    
      %ecef COORDINATE for sat 1
      x_ECEFcor(k)=coordinate_ECEF(1);
      y_ECEFcor(k)=coordinate_ECEF(2);
      z_ECEFcor(k)=coordinate_ECEF(3);
      x_cor1(k)=Xeci1;
          y_cor1(k)=Yeci1;
          z_cor1(k)=Zeci1;
          %ecef COORDINATE for sat 2
      x_ECEFcor1(k)=coordinate_ECEF1(1);
      y_ECEFcor1(k)=coordinate_ECEF1(2);
      z_ECEFcor1(k)=coordinate_ECEF1(3);
          x_cor2(k)=Xeci2;
          y_cor2(k)=Yeci2;
          z_cor2(k)=Zeci2;
          %ecef COORDINATE for sat 3
      x_ECEFcor2(k)=coordinate_ECEF2(1);
      y_ECEFcor2(k)=coordinate_ECEF2(2);
      z_ECEFcor2(k)=coordinate_ECEF2(3);
          x_cor3(k)=Xeci3;
          y_cor3(k)=Yeci3;
          z_cor3(k)=Zeci3;
          %ecef COORDINATE for sat 4
      x_ECEFcor3(k)=coordinate_ECEF3(1);
      y_ECEFcor3(k)=coordinate_ECEF3(2);
      z_ECEFcor3(k)=coordinate_ECEF3(3);
%% this commented part plot the satellite motion in ECEF  
%  ECEFarray(k-1)=plot3 (x_ECEFcor(k), y_ECEFcor(k), z_ECEFcor(k),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
%    hold on
% %      ECEFarray1(k-1)=plot3 (x_ECEFcor1(k), y_ECEFcor1(k), z_ECEFcor1(k),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','g','MarkerSize', 6);
%    hold on
% %      ECEFarray2(k-1)=plot3 (x_ECEFcor2(k), y_ECEFcor2(k), z_ECEFcor2(k),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
%     hold on
% %      ECEFarray3(k-1)=plot3 (x_ECEFcor3(k), y_ECEFcor3(k), z_ECEFcor3(k),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 6);
%     hold on     
%% this part plot four satellite motion in ECI frame
if k~=1   
    array=plot3 (x_cor(k-1), y_cor(k-1), z_cor(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 6);
    hold on
    array1=plot3 (x_cor1(k-1), y_cor1(k-1), z_cor1(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','g','MarkerSize', 6);
    hold on
    array2=plot3 (x_cor2(k-1), y_cor2(k-1), z_cor2(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','b','MarkerSize', 6);
    hold on
%        position(k-1)=line([0 x_cor(k-1)],[0 y_cor(k-1)], [0 z_cor(k-1)],'Color', 'black', 'LineWidth', 0.2);%
     array3=plot3 (x_cor3(k-1), y_cor3(k-1), z_cor3(k-1),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 6);
%plot label 
    legend([array,array1,array2,array3],'Sat 1','Sat 2','Sat 3','Sat 4')     
end
%% this commented part could plot radial vector of satellite not necessary depends on user demand
%    position1(k)=line([0 x_cor1(k-1)],[0 y_cor1(k-1)], [0 y_cor1(k-1)],'Color', 'black', 'LineWidth', 0.2);%    
%     position2(k)=line([0 x_cor2(k-1)],[0 y_cor2(k-1)], [0 z_cor2(k-1)],'Color', 'black', 'LineWidth', 0.2);%    
%      position3(k)=line([0 x_cor3(k-1)],[0 y_cor3(k-1)], [0 z_cor2(k-1)],'Color', 'black', 'LineWidth', 0.2);%    
%visualize satellite movement on the 2D earth surface
%calculate latitude and longtitude coordinates
    lat(k-1)=atan(z_cor(k-1)/sqrt(x_cor(k-1)*x_cor(k-1)+y_cor(k-1)*y_cor(k-1)))*(180/pi);
    lat1(k-1)=atan(z_cor1(k-1)/sqrt(x_cor1(k-1)*x_cor1(k-1)+y_cor1(k-1)*y_cor1(k-1)))*(180/pi);
    lat2(k-1)=atan(z_cor2(k-1)/sqrt(x_cor2(k-1)*x_cor2(k-1)+y_cor2(k-1)*y_cor2(k-1)))*(180/pi);
    lat3(k-1)=atan(z_cor3(k-1)/sqrt(x_cor3(k-1)*x_cor3(k-1)+y_cor3(k-1)*y_cor3(k-1)))*(180/pi);
    ECIX=[cos(GMST) sin(GMST) 0];
    Pos=[real(x_cor(k-1)) real(y_cor(k-1)) 0];
    cvec = cross(ECIX,Pos);
    Pos1=[real(x_cor1(k-1)) real(y_cor1(k-1)) 0];
    cvec1 = cross(ECIX,Pos1);
    Pos2=[real(x_cor2(k-1)) real(y_cor2(k-1)) 0];
    cvec2 = cross(ECIX,Pos2);
    Pos3=[real(x_cor3(k-1)) real(y_cor3(k-1)) 0];
    cvec3 = cross(ECIX,Pos3);
    angleyz1 = mod(sign(dot([0 0 1],cvec1))*atan2(norm(cvec1),dot(ECIX,Pos1)),2*pi);
    long1(k-1) =(180/pi)* angleyz1;
     angleyz3 = mod(sign(dot([0 0 1],cvec3))*atan2(norm(cvec3),dot(ECIX,Pos3)),2*pi);
    long3(k-1) =(180/pi)* angleyz3;
    angleyz2 = mod(sign(dot([0 0 1],cvec2))*atan2(norm(cvec2),dot(ECIX,Pos2)),2*pi);
    long2(k-1) =(180/pi)* angleyz2;
    angleyz = mod(sign(dot([0 0 1],cvec))*atan2(norm(cvec),dot(ECIX,Pos)),2*pi);
    long(k-1) =(180/pi)* angleyz;
        %calculate geodetic coordinate from ECEF coordinates
        long_ecef(k-1)=atan2(real(coordinate_ECEF(2)),real(coordinate_ECEF(1)))*180/pi;
        if long_ecef(k-1)<0
            long_ecef(k-1)=360-abs(long_ecef(k-1))
        end
        lat_ecef(k-1)=asin(real(coordinate_ECEF(3))/(r_magnitude(k-1)-EARTH_RAIDUS*e^2))*180/pi;
%         if lat_ecef(k-1)<0
%             lat_ecef(k-1)=360-abs(lat_ecef(k-1))
%         end
        long_ecef1(k-1)=atan2(real(coordinate_ECEF1(2)),real(coordinate_ECEF1(1)))*180/pi;
        if long_ecef1(k-1)<0
            long_ecef1(k-1)=360-abs(long_ecef1(k-1));
        end
        lat_ecef1(k-1)=asin(real(coordinate_ECEF1(3))/(r1_magnitude(k-1)-EARTH_RAIDUS*e^2))*180/pi;
        
        long_ecef2(k-1)=atan2(real(coordinate_ECEF2(2)),real(coordinate_ECEF2(1)))*180/pi;
        if long_ecef2(k-1)<0
            long_ecef2(k-1)=360-abs(long_ecef2(k-1));
        end
        lat_ecef2(k-1)=asin(real(coordinate_ECEF2(3))/(r2_magnitude(k-1)-EARTH_RAIDUS*e^2))*180/pi;
        
        long_ecef3(k-1)=atan2(real(coordinate_ECEF3(2)),real(coordinate_ECEF3(1)))*180/pi;
        if long_ecef3(k-1)<0
            long_ecef3(k-1)=360-abs(long_ecef3(k-1))
        end
        lat_ecef3(k-1)=asin(real(coordinate_ECEF3(3))/(r3_magnitude(k-1)-EARTH_RAIDUS*e^2))*180/pi;
%         if lat_ecef3(k-1)<0
%             lat_ecef3(k-1)=360-abs(lat_ecef3(k-1))
%         end 
    end
    if k~=2
%     line([x_cor(k-1) x_cor(k)],[y_cor(k-1) y_cor(k)],[z_cor(k-1) z_cor(k)],'Color', 'red', 'LineWidth', 2);
%     line([x_cor1(k-1) x_cor1(k)],[y_cor1(k-1) y_cor1(k)],[z_cor1(k-1) z_cor1(k)],'Color', 'green', 'LineWidth', 2);
%     line([x_cor2(k-1) x_cor2(k)],[y_cor2(k-1) y_cor2(k)],[z_cor2(k-1) z_cor2(k)],'Color', 'blue', 'LineWidth', 2);
%     line([x_cor3(k-1) x_cor3(k)],[y_cor3(k-1) y_cor3(k)],[z_cor3(k-1) z_cor3(k)],'Color', 'yellow', 'LineWidth', 2);
   
   
    end
   %update mean anomaly
   Mean_anomaly=Mean_anomaly+sqrt((GM/(a*a*a)))*satrunspeed;
end
%plot ECEF projection 
figure (5)
set(gcf,'Menubar','none','Name','Satellite Projection Earth 2D surface ', ... 
    'NumberTitle','off','Position',[10,350,1000,500], ... 
    'Color',[1 1 1]); 
hold on
image([0 360],[-90 90],topo,'CDataMapping', 'scaled');
colormap(topomap1);
axis equal
axis ([0 360 -90 90]);
plot (167.717,8.717,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
plot (360-76.496, 42.440,'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 10);
for i=1:k-3
 %plot the projection on earth 2D surface in ECEF frame
     a=plot (long_ecef(i),lat_ecef(i),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 7);
   % plot (long_ecef1(i),lat_ecef1(i),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 4);
     b=plot (long_ecef2(i),lat_ecef2(i),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','y','MarkerSize', 7);
    % plot (long_ecef3(i),lat_ecef3(i),'s', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 4);
    legend([a,b],'Sat 1','Sat 2')
    xlabel('Latitude')
    ylabel('Longtitude')
     pause (0.1);
end
 
 hold off
 title('Communication Satellite Projection')
 %% plot satellite orbital height variation
 figure(4)
 x=linspace(0,k-2,k-1);
 x=x*orbital_period/length(x);
 plot(x(1:k-1),satelliteheight(1:k-1),'LineWidth',2)
 hold on
 plot(x(1:k-1),satelliteheight1(1:k-1),'LineWidth',2)
 hold on
 plot(x(1:k-1),satelliteheight2(1:k-1),'LineWidth',2)
 hold on
 plot(x(1:k-1),satelliteheight3(1:k-1),'LineWidth',2)
 hold on
 grid on
 legend('Satellite 1 Orbit Height','Satellite 2 Orbit Height','Satellite 3 Orbit Height','Satellite 4 Orbit Height')
 ylabel('Orbital Height(unit:m)')
end
% 0.5*v_minimum^2=GM/earth_radius^2;



          