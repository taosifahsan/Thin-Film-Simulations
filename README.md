# Thin-Film-Simulations
We are simulating retraction of a thin film after rupture in various situations. 

There are three source codes. All are written in matlab. They are implemented from scratch without any use of library or packages. 

# Source Code 1: thin_film.m 
We are simulating evolution of a finite thin film in laboratory frame. Section III used this.

  The following parameters controls the motion of the thin film.
  
  Oh = Ohnesorge number (Oh = viscosity/sqrt(density*surface tension co-efficient*initial height).
  
  epsilon = aspect ratio of the thin film (initial length/initial height).
  
  N = number of bins in x-axis.
  
  T = total time of simulation.
  
  NT = number of bins in time-axis.
  
  Nframe = number of frames to be saved.
  
  initial_condition(epsilon,N,1,0) initializes the thin film.
  
  simulate(x_0,v_0,h_0,Oh,T,NT,Nframe) simulates the thin film.
  
After getting the results, plot using the following methods. 
  
  plot_heights(t,x,h,epsilon,Oh)    this plots height profile.
  
  plot_velocities(t,x,v,epsilon,Oh)  this plots velocity profile.
  
  plot_tip_velocity(t,v,x,h,epsilon,Oh) this plots retraction speed of the tip.
  
  plot_maximum_height(t,h,epsilon,Oh)  this plots the maximum height in the thin film.
  
# Source Code 2: thin_film_L=infty.m

We are simulating evolution of an infinite thin film in tip frame. Section VI used this.
  
  The following parameters controls the motion of the thin film.
  
  Oh = Ohnesorge number (Oh = viscosity/sqrt(density*surface tension co-efficient*initial height).
  
  L0 = the part of x-axis we are simulating.
  
  N = number of bins in x-axis.
  
  T = total time of simulation.
  
  NT = number of bins in time-axis.
  
  Nframe = number of frames to be saved.
  
  [t,x,v,h]=simulate_thin_film(Oh,h0,L0,N,T,NT) simulates the thin film.
  
After getting the results, plot using the following methods.
  
  plot_heights(t,x,h,Oh)  this plots height profile.
  
  plot_velocities(t,x,v,Oh) this plots velocity profile.
  
  plot_tip_velocity(t,v,Oh) this plots retraction speed of the tip.
  
  plot_maximum_height(t,h,Oh) this plots the maximum height in the thin film.
  
# Source Code 3: thin_film_Oh>>1.m.

We used this for numerical calculation of a more general finite thin film with Oh>>1 for various aspect-ratio/Oh values. Subsection V.7 used this. 

  L = aspect-ratio/Oh;
  
  T = total time of simulation;
  
  N = number of spatial bins;
  
  NT = numbers of time bins;

Simulate and plot the results for an array of parameters using the following function. Setting plot_all = 1 will plot all of the results of the simulation. Otherwise it will only plot the v_max vs L, t_max vs L, t_phase vs L and t_phase_duration vs L graphs. 

  run_code_array(L,T,N,NT,plot_all)

# Source Code 3: thin_film_L=infty_Oh>>1.m

We used this for numerical calculation of infinite thin film with Oh>>1 and aspect-ratio/Oh>>1. Subsection V.7 used this. 

  N =  number of space bins.  
  
  L = region of x-axis we are simulating.
  
  dx = L/N bin gap.
  
  x = linspace(0,L,N) x-axis.
  
  NT = number of time bins.
  
  T = total time.
  
  dt = T/NT time gap.

  Nframe = number of frames to be saved.

  [t,h,u] = simulate(T,dt,dx,L,Nframe) this calculates the results.

  Ts = round(linspace(1,5,5)*length(t_r)/5) plot these snapshots during profile plot.
  
  fsize = 18 fontsize.
  
Plot the results using these methods.
  
  plot_heights(t_r,x,h,Ts,fsize).
  
  plot_velocities(t_r,x,h,Ts,fsize). 
  
  plot_V_TC(t_r,u,fsize). 


# Extra Source Codes: 

There are some additional programs. We calculated traffic flow problem using traffic_flow.m in appendix A.

We compared the runtime of thomas algorithm and direct matlab matrix inversion using tridiagonal_alg_check.m in appendix B.


  
