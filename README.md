
# Project: **Orbit simulation for formation flying in LEO**.<br/>
Matlab codes and thesis for Master project at the university of Bremen, department of Communications engineering. <br/> 
<br/>
The basic introduction and partial research results and demonstrations are given here. <br/>
<br/>
Main research goals in this project are summarized as follows:<br/>
           1. Design a satellite orbit propagator in low earth orbit to simulate satellite motion with assistance of orbital elements <br/>
           2. Satellite motion simulation and transform on different frames. ECI frame, ECEF and LVLH frame.<br/>
           3. Satellite relative motion dynamics with HCW equations<br/>
           4. Kalman filter implementation based on HCW equations and sensor measurements for satellite position and velocity estimation to achieve satellite naviagation and optimal control<br/>
           5. The extensive investigation about Drag-modified HCW equations and J2-modified HCW equations, which involves perturbation of atmospheric drag force and earth flattening effect.<br/>
           6. Extended Kalman filter implementation based on these modified HCW equation and sensor measurements for satellite position and velocity estimation  under above mentioned perturbations<br/>

<br/>
The illustration of a Deputy- Chief satellite network on LVLH(Local vertical local horizontal) frame  

<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900189-8cae08c1-619d-460c-87fe-bf5344676705.png"
" width="400" />
  </p> 
The HCW equation for inter-satellite motion dynamics:
<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900413-94bc1ef8-395e-47c6-b3e6-27babf01fdcb.png"
" width="400" />
  </p> 



 Algorithm structure of Kalman filter to recursively estimate state of dynamic system for linear system.            
<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900943-f9cbc68c-d451-4304-90a6-686ff666f8a1.png"
 " width="400" />
</p> 
System dynamics model and discrete system state update which are derived from HCW equation for position and velocity estimation on LVLH frame for deputy satellite
<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900944-2bcf9db7-d48b-49b0-aede-b15fe2409aaf.png"
" width="400" />
  </p> 
<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900946-a675ed79-784e-4c9e-93a4-09d20aeac70a.png"
" width="400" />
  </p> 
Kalman filter for estimation of velocity and position of deputy satellite with respect to chief satellite. The estimation performance is evaluated with different level of measurement noise variance ($ \Sigma_w^2$). 

<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900780-3ece325b-9838-4714-9403-0632797d8a23.png"
" width="800" />
  </p> 
This sentence uses `$` delimiters to show math inline:  $\sqrt{3x-1}+(1+x)^2$
