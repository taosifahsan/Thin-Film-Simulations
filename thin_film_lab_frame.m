%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Main
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% set parameters
%--------------------------------------------------------------------------
Oh = 10; %Ohnesorge number (Oh = μ/sqrt(ργh_0))
epsilon = 10; %aspect ratio of the thin film
N = 2000; % number of space bins
T = 0.5*(4*Oh); % total time in v_t + v*v_x = 4 Oh (hv_x)_x + k_x equation
NT = round(10^5);% number of time bins
Nframe = min(NT,10^4); % number of frames to be saved
%--------------------------------------------------------------------------
% simulate
%--------------------------------------------------------------------------
[x_0,v_0,h_0]=initial_condition(epsilon,N,1,0); %initialize
[t,x,v,h] = simulate(x_0,v_0,h_0,Oh,T,NT,Nframe);%simulate

%--------------------------------------------------------------------------
% plot the simulations
%--------------------------------------------------------------------------
%height and velocity profile plotted
plot_heights(t,x,h,epsilon,Oh)    
plot_velocities(t,x,v,epsilon,Oh) 
plot_tip_velocity(t,v,x,h,epsilon,Oh)
plot_maximum_height(t,h,epsilon,Oh)     


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Define Initial Conditions Here
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Intialize position, velocity and height
% [x,v,h] = initial_condition(aspect ratio, number of time bins,
% height scale,velocity scale)
%--------------------------------------------------------------------------
function [x,v,h] = initial_condition(epsilon,N,h0,v0)
    x = linspace(0,epsilon*h0,N); %initialize positions
    v = v0*linspace(0,0,N); %initialize velocity
    L0 = h0*epsilon;
    xp = L0-x;
    a = 1/20;
    f = 1/2-a-xp+1/2*sqrt((1+2*a)^2+4*xp.*(xp+2*a-1));
    h = h0*sqrt(1-f.^2); %initialize heights
    %h = h0*tanh(sqrt(epsilon-x));
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Simulate The Thin Film using Crank Nicolson and Picard Iteration Method
%--------------------------------------------------------------------------
% first simulated using v_t + v*v_x = 4 Oh (hv_x)_x + k_x equation
% then converted back to V_T + V*V_x = (HV_X)_X + K_X using 
% T = t/(4Oh), X = x/(4Oh)

% [time, position, velocity, height]=simulate(initial position,
% velocity, initial height, Ohnesorge Number, Total Time, 
% Number of Time Bins, Number of space bins, Number of saved frames);

% total_error ~ O(dx^2+dt^2) where dx, dt is spatial and temporal step size 
% runtime ~ O(N*NT) % N, NT number of spatial and temporal bins
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [t,x,v,h] = simulate(x_0,v_0,h_0,Oh,T,NT,Nframe)    
    %the inital position, velocity and height of space bins
    x_n=x_0;
    v_n=v_0;
    h_n=h_0;

    N = length(x_0); %number of position bins

    N_write_time = 100; %100 checkpoints in the terminal for sanity check
    elapsed_time = 0; %run time starts
    tic 
    
    dt = T/NT; %time step size

    % the information will be saved in these matrices
    % we only save a fraction of the time data to save space
    % N_frame (No of saved time bins) < NT (No of time bins)
    frame_step = round(NT/Nframe); % time bins between two saved frame
    h = zeros(N,Nframe);
    t = linspace(0,0,Nframe);
    v = zeros(N,Nframe);
    x = zeros(N,Nframe);
    
    i_check = 0; % variable related check point to be printed on terminal
    Total_repeat = 0;
    for it = 1:NT 
        % rough explicit prediction for future velocity, position and
        % height
        v_f = v_step_explicit(h_n,v_n,x_n,dt,Oh);
        x_f = x_n+(v_n+v_f)/2*dt;
        h_f = h_step(h_n,v_n,x_n,v_f,x_f,dt);

        % crank nicolson prediction for height using explicitly predicted
        % future velocity using picard iteration
        MSE = 1e+12;  % initial magnitude squared error 
        convergence_cutoff = 1e-9; %required cutoff for valid solution
        Num_repeats = 0;
        Num_Failure = 20;
        % iterate until convergence achieved or convegence fails 
        % within certain steps 
        while(MSE>convergence_cutoff && Num_repeats<Num_Failure) 
            v_temp = v_f;
            % a better prediction of v using crank nicolson method
            v_f = v_step_crank_nicolson(h_n,v_n,x_n,h_f,x_f,dt,Oh);
            MSE = sum((v_temp-v_f).^2); % calculate MSE between iterations
            
            % prediction of position height
            x_f = x_n+(v_n+v_f)/2*dt; 
            h_f = h_step(h_n,v_n,x_n,v_f,x_f,dt);

            %increase number of repeats
            Num_repeats = Num_repeats+1;
        end
        Total_repeat = Total_repeat + Num_repeats;
        
        % assert convergence, break loop if convergence fails
        if (Num_repeats==Num_Failure) 
            fprintf("CONVERGENCE FAILED, MSE = %5.4e\n = ", MSE);
            break
        end


        % save everything here with reduced resolution to save space
        if mod(it,frame_step) == 0
            t(round(it/frame_step))=it/NT*T;
            h(:,round(it/frame_step))=h_n;
            x(:,round(it/frame_step))=x_n;
            v(:,round(it/frame_step))=v_n;
        end
        
        % future data becomes present data
        x_n = x_f;
        v_n = v_f;
        h_n = h_f;
        
        % check points during runtime, print checks in terminal
        if mod(it,round(NT/N_write_time))==0
            elapsed_time_past = elapsed_time;
            elapsed_time = toc;
            total_time = elapsed_time-elapsed_time_past;
            time = it/NT*T;
            i_check = i_check+1;

            sentence = ['%d. t = %2.4f, L = %5.4f, v(L)=%5.4f,' ...
                '  h(0)= %2.4f, sim = %2.2f, run time = %2.1fs,' ...
                ' Avg Repeat = %2.2f\n'];
            fprintf(sentence, i_check, time/(4*Oh), x_n(N), -v_n(N), h_n(1), ...
                x_n(N)*h_n(1)/x_0(N), total_time, Total_repeat/(NT/N_write_time))

            Total_repeat = 0;
        end


    end
    toc
    t = t/(4*Oh);
    x = x/(4*Oh);
end

%--------------------------------------------------------------------------
% Predict future v using explicit method
%--------------------------------------------------------------------------
function v_f = v_step_explicit(h_n,v_n,x_n,dt,Oh)
    hx_n = D(h_n,x_n);
    kx_n = kappa_x(h_n,hx_n,x_n);
    
    vx_n = D(v_n,x_n);
    vxx_n = D2(v_n,x_n);
   
    %future velcity
    dv_dt = 4*Oh*(hx_n.*vx_n./h_n+vxx_n)+kx_n;
    v_f = v_n+dt*dv_dt; 

    %boundary conditions
    v_f(1)=0; % v(0)=0
    v_f(length(v_f))=v_f(length(v_f)-1); %v_x(L)=0
end

%--------------------------------------------------------------------------
% Predict Future v using crank nicolson method
%--------------------------------------------------------------------------
function v_f = v_step_crank_nicolson(h_n,v_n,x_n,h_f,x_f,dt,Oh)
    N = length(x_n);
    
    % height derivatives of present
    hx_n = D(h_n,x_n);

    % curvature derivative of present
    kx_n = kappa_x(h_n,hx_n,x_n);
    
    % contribution of explicit part on (hv_x)_x/h+k_x
    vx_n = D(v_n,x_n);
    vxx_n = D2(v_n,x_n);
    vt_exp = 4*Oh*(hx_n.*vx_n./h_n+vxx_n)+kx_n;

    %create matrix and vector for v for PDE equations needed for implicit
    hx_f = D(h_f,x_f);

    % create matrix A for implicit part
    % these three digonal arrays are only non zero entries of Matrix A
    A_0 = linspace(0,0,N); %A(i,i) array
    A_plus = linspace(0,0,N); %A(i,i+1) array
    A_minus = linspace(0,0,N); %A(i,i-1) array

    for ix = 2:N-1
        dx_f = (x_f(ix+1)-x_f(ix-1))/2;
        fact_A = 4*Oh*dt/dx_f^2;
        fact_B =hx_f(ix)/h_f(ix)*dx_f/2;
        A_0(ix)=1+fact_A ; %update the diagonal entries

        %update two other off diagonal entries
        A_plus(ix) =-(1+fact_B)*fact_A/2; 
        A_minus(ix)=-(1-fact_B)*fact_A/2;
    end     
   
    %enforce boundary conditions for v on the Matrix
    %first boundary condition v(0)=0
    A_0(1)=1; 
    %Second boundary condition v_x(L) = 0
    A_minus(N)=-1;
    A_0(N)=1;
    
    %curvature derivative in future   
    kx_f = kappa_x(h_f,hx_f,x_f);
        
    % crank nicolson vector
    d = v_n+dt*kx_f/2+dt*vt_exp/2;
        
    %enforce boundary conditions of v on the vector
    %first boundary condition v(0)=0
    d(1)=0; 
    %second boundary condition v_x(L) = 0
    d(N)=0; 
 
    % future velocity prediction
    v_f=solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

%--------------------------------------------------------------------------
% Predict Future h using crank nicolson method
%--------------------------------------------------------------------------
function h_f = h_step(h_n,v_n,x_n,v_f,x_f,dt)
    N = length(x_n);
    h_f = h_n;
    for ix = 2:N-2
            % central finite difference of position upto third order
            % finite difference of velocity upto two order
            dvi = (v_n(ix+1)-v_n(ix-1))/2;
            dxi = (x_n(ix+1)-x_n(ix-1))/2;
            dxf = (x_f(ix+1)-x_f(ix-1))/2;
            dvf = (v_f(ix+1)-v_f(ix-1))/2;
            h_f(ix)=h_n(ix)/(1+dvf/dxf*dt/2)*(1-dvi/dxi*dt/2);
    end  
    
    %boundary conditions
    h_f(N)=0; %h(L) = 0
    
    ratio = sqrt((x_f(N)-x_f(N-1))/(x_f(N)-x_f(N-2)));
    h_f(N-1)=h_f(N-2)*ratio; %lim x->L h(x) ~ (L-x)^(1/2)

    h_f(1)=h_f(2); %h_x(0) = 0
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Extra Methods for The Algorithm
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Calculate curvature derivative function
%--------------------------------------------------------------------------
function kx = kappa_x(h,hx,x)
    hxx = D2(h,x);
    k = hxx./(1+hx.^2).^(3/2); % curvature
    kx = D(k,x); % curvature derivative
end

%--------------------------------------------------------------------------
% Tridiagonal Thomas Algorithm for Matrix Inversion
%--------------------------------------------------------------------------
function x = solve_tridiagnoal(A_0,A_minus,A_plus,d)
    n = length(d);
    x = linspace(0,0,n);
   
    for i = 2:n
        w = A_minus(i)/A_0(i-1);
        A_0(i)=A_0(i)-w*A_plus(i-1);
        d(i)=d(i)-w*d(i-1);
    end

    x(n)=d(n)/A_0(n);
    for i = n-1:-1:1
        x(i)=(d(i)-A_plus(i)*x(i+1))/A_0(i);
    end

end

%--------------------------------------------------------------------------
% First partial derivative wrt position, finite central difference 
%--------------------------------------------------------------------------
function y=D(f,x)
    N=length(f);
    df=linspace(0,0,N);
    dx=linspace(1,1,N);
    for i = 2:N-1
        df(i)=(f(i+1)-f(i-1))/2;
        dx(i)=(x(i+1)-x(i-1))/2;
    end
    y=df./dx;
    y(1)=2*y(2)-y(3);
    y(N)=2*y(N-1)-y(N-2);
end 

%--------------------------------------------------------------------------
% Second partial derivative wrt position, finite central difference  
%--------------------------------------------------------------------------
function y=D2(f,x)
    N=length(f);
    y=linspace(0,0,N);
    for i = 2:N-1
        dx = (x(i+1)-x(i-1))/2;
        d2f=f(i+1)+f(i-1)-2*f(i);
        y(i)=d2f/dx^2;
    end
    y(1)=2*y(2)-y(3);
    y(N)=2*y(N-1)-y(N-2);
end 


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Plotting Methods
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Plot Heights
%--------------------------------------------------------------------------
function plot_heights(t,x,h,epsilon,Oh)  
    cmap = colororder(); %color sequence of time snaps
    Nframe = length(t);
    Ts = round([1,linspace(1,4,4)*Nframe/4]); % timesnaps showed in figure 
    figure 
    for i = 1:length(Ts)
        label_name = 'T = '+string(round(t(Ts(i)),2));
        if i ==1
            label_name = sprintf('T = %2.1e',t(Ts(1)));
        end 
        plot(x(:,Ts(i)),h(:,Ts(i)),'color', cmap(i,:), ...
            'linewidth',2, DisplayName=label_name);
        hold on
        plot(x(:,Ts(i)),-h(:,Ts(i)),'color', cmap(i,:),'linewidth',2, ...
            'HandleVisibility','off');
        hold on
    end 
    label_figure('Position, X', 'Height, H',epsilon,Oh);
    hold off
end


%--------------------------------------------------------------------------
% Plot Velocities
%--------------------------------------------------------------------------
function plot_velocities(t,x,v,epsilon,Oh) 
    NT = length(t);
    Ts = round([1,linspace(1,4,4)*NT/4]); % timesnaps showed in figure 
    figure 
    for i = 1:length(Ts)
        label_name = 'T = '+string(round(t(Ts(i)),2));
        if i ==1
            label_name = sprintf('T = %2.1e',t(Ts(1)));
        end 
        
        plot(x(:,Ts(i)),-v(:,Ts(i)),'linewidth',2, ...
            DisplayName=label_name);
        hold on
    end     
    label_figure('Position, X', 'Velocity, V(X,T)',epsilon,Oh);
    hold off
end

%--------------------------------------------------------------------------
% Plot Velocity of Free End Point
%--------------------------------------------------------------------------
function plot_tip_velocity(t,v,x,h,epsilon,Oh)
    N = length(v(:,1));
    vmax = linspace(0,0,length(t));
    for it = 1:length(t)
        hmax = h(1,it);
        ix = round(N*(x(N,it)-hmax/(4*Oh))/x(N,it));
        v_cm = sum(v(ix:1:N,it).*h(ix:1:N,it))/sum(h(ix:1:N,it));
        vmax(it)=v_cm;
    end
    figure
    plot(t,-vmax,'linewidth', 2, DisplayName='Numerical V_{tip}'); 
    label_figure('time, T', 'Tip Velocity, U(T)',epsilon,Oh);
    hold off
end


%--------------------------------------------------------------------------
% Plot Maximum Height
%--------------------------------------------------------------------------
function plot_maximum_height(t,h,epsilon,Oh)     
    hmax = linspace(0,0,length(t));
    for it=1:length(t)
        hmax(it)=max(h(:,it));
    end
    figure
    plot(t,hmax,'linewidth', 2, DisplayName='Numerical'); 
    label_figure('time, T ', 'Maximum Height, H_{max}(T)',epsilon,Oh);
end

%--------------------------------------------------------------------------
% Labeling the Figures
%--------------------------------------------------------------------------
function label_figure(x_label,y_label,epsilon,Oh)
    % title of the calculations and figures
    fontsize = 15; %fontsize 
    title_string = 'Oh = '+string(Oh)+', $\varepsilon$ = '+string(epsilon);

    xlabel(x_label,'fontsize',fontsize,'fontname','times') %x-axis label
    ylabel(y_label,'fontsize',fontsize,'fontname','times') %y-axis label

    title(title_string,'fontsize',fontsize,'Interpreter','latex') %titile
    legend('Location','best','Fontsize',fontsize,'fontname','times')

    grid on
    grid minor
    set(gca,'FontSize',fontsize,'fontname','times')
end
