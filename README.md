
# Project: **Orbit simulation for formation flying in LEO**.<br/>
Matlab codes for Master project at the university of Bremen, department of Communications engineering. <br/>

Main contributions in this project are summarized as follows:<br/>
           1. Design a satellite orbit propagator in low earth orbit to simulate satellite motion with assistance of orbital elements <br/>
           2. Satellite motion simulation and transform on different frames. ECI frame, ECEF and LVLH frame.<br/>
           3. Satellite relative motion dynamics with HCW equations<br/>
           4. Kalman filter implementation based on HCW equations and sensor measurements for satellite position and velocity estimation to achieve satellite naviagation and optimal control<br/>
           5. The extensive investigation about Drag-modified HCW equations and J2-modified HCW equations, which involves perturbation of atmospheric drag force and earth flattening effect.<br/>
           6. Extended Kalman filter implementation based on these modified HCW equation and sensor measurements for satellite position and velocity estimation  under above mentioned perturbations<br/>

<br/>
The illustration of a Deputy Chief satellite network on LVLH(Local vertical local horizontal) frame  

<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900189-8cae08c1-619d-460c-87fe-bf5344676705.png"
" width="400" />
  </p> 
The HCW equation for inter-satellite motion dynamics:
<p align="center">
  <img src="https://user-images.githubusercontent.com/89796179/198900413-94bc1ef8-395e-47c6-b3e6-27babf01fdcb.png"
" width="400" />
  </p> 
