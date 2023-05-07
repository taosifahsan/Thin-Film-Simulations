Oh = 0.1; %Ohnesorge number (Oh = μ/sqrt(ργh_0))
L0 = 10
h0 = 1;
%here height or h is half width of the thin film

dx = L0/N; %bin gap

T = 20; %total time
dt = T/NT; %time gap

Nframe = min(NT,10000); %number of frames to be saved

% [time, position, velocity, height]=
% simulate_thin_film(Ohnesorge Number, ratio of initial length and height,
% Total Time, time steps,position steps,nature of position distribution,
% initial length,number of frames of saved data);

[t,x,v,h]=simulate_thin_film(Oh,h0,L0,N,T,NT);
plot_heights(t,x,h,Oh)  
plot_velocities(t,x,v,Oh) 
plot_tip_velocity(t,v,Oh)
plot_maximum_height(t,h,Oh)     




%--------------------------------------------------------------------------
% Intial height
%--------------------------------------------------------------------------
function y = h_in(x,h0,Oh)
    %x = x/(4*Oh);
    %a = 1/100;
    %y = 1/2-a-x+1/2*sqrt((1+2*a)^2+4*x.*(x+2*a-1));
    %y = h0*sqrt(1-y.^2);

    y = linspace(1,1,length(x));
    y(1)=0;
    y(2)=sqrt(1/2);


    %y=tanh((x/h0).^(1/2))*h0;
end

%--------------------------------------------------------------------------
% Simulate The Thin Film using Crank Nicolson Method
%--------------------------------------------------------------------------


function [t,x,v,h] = simulate_thin_film(Oh,h0,L0,N,T,NT)
    %the intial position of space bins

    x = linspace(0,L0,N);
    dx = x(N)/N;
    dt = T/NT;

    %the inital velocity and height of space bins
    v_n = zeros(N,1)'; 
    h_n = h_in(x,h0,Oh);

    N_write_time = 100;
    elapsed_time = 0;
    tic 
    % the information will be saved here
    % we only save a fraction of the time data to save space
    % N_frame (No of saved time bins) < N_time (No of time bins)

    Nframe = min(NT,10000); %number of frames to be saved
    frame_step = round(NT/Nframe); % time bins between two saved frame
    h = zeros(N,Nframe);
    t = linspace(0,0,Nframe);
    v = zeros(N,Nframe);
    
    i_check = 0; % variable related check point to be printed on terminal
    Total_repeat = 0;
    for it = 1:NT
        % rough explicit prediction for future velocity and position
        v_f =  v_step_explicit(h_n,v_n,dx,dt,Oh);
        h_f =  h_step(h_n,v_n,v_f,dx,dt);
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

            v_f = v_step_crank_nicolson(h_n,v_n,h_f,v_f,dx,dt,Oh);
            MSE = sum((v_temp-v_f).^2); % calculate MSE between iterations

            h_f = h_step(h_n,v_n,v_f,dx,dt);
            Num_repeats = Num_repeats+1;
        end
        Total_repeat = Total_repeat + Num_repeats;
        % assert convergence, break loop if convergence fails
        if (Num_repeats==Num_Failure) 
            fprintf("CONVERGENCE FAILED, MSE = %5.4e\n = ", MSE);
            break
        end


        % save everything here
        if mod(it,frame_step) == 0
            t(round(it/frame_step))=it/NT*T;
            h(:,round(it/frame_step))=h_n;
            v(:,round(it/frame_step))=v_n;
        end
        
        % future data becomes present data
        v_n = v_f;
        h_n = h_f;
        
        %check point
        if mod(it,round(NT/N_write_time))==0
            elapsed_time_past = elapsed_time;
            elapsed_time = toc;
            total_time = elapsed_time-elapsed_time_past;
            time = it/NT*T;
            i_check = i_check+1;
            hmax = max(h_n);
            hx_N = (h_n(N)-h_n(N-1))/dx;
            sentence = ['%d. t = %2.4f, v_tip=%5.4f, hmax = %5.4f,' ...
                ' hx_wall = %5.4f, run time = %2.1fs, Avg Repeat = %2.2f\n'];
            fprintf(sentence, i_check, time, v_n(1), hmax, hx_N, total_time, ...
                Total_repeat/(NT/N_write_time))
            Total_repeat = 0;
        end
    end
    toc
end

%--------------------------------------------------------------------------
% Predict Future v
%--------------------------------------------------------------------------
function v_f_cnm = v_step_crank_nicolson(h_n,v_n,h_f,v_f,dx,dt,Oh)
    N = length(h_n);
    
    % height derivatives of present
    hx_n = D(h_n,dx);
    % height derivatives of future
    hx_f = D(h_f,dx);

    % curvature derivative of present
    k_n = kappa(h_n,hx_n,dx,Oh);
    kx_n = D(k_n,dx); % curvature derivative
    % curvature derivative in future   
    k_f = kappa(h_f,hx_f,dx,Oh);
    kx_f =D(k_f,dx);
    
    % contribution of explicit part on (hv_x)_x/h+k_x
    vx_n = D(v_n,dx);
    vxx_n = D2(v_n,dx);

    %factor in front of vx
    A_n = hx_n./h_n+(v_n(1)-v_n); %present
    A_f = hx_f./h_f+(v_f(1)-v_f); %future

    % crank nicolson vector
    d = v_n+dt/2*(kx_f+kx_n+A_n.*vx_n+vxx_n);
        
    %enforce boundary conditions of v on the vector
    %first boundary condition v(0)=0
    
    %second boundary condition v_x(L) = 0
    

    % create matrix A for implicit part
    R = dt/dx^2;
    % these three digonal arrays are only non zero entries of Matrix A
    
    A_0 = zeros(N,1)'+(1+R);
    A_plus =-R/2*(1+A_f/2*dx);
    A_minus=-R/2*(1-A_f/2*dx); 
   
    %enforce boundary conditions for v on the Matrix
    %first boundary condition v_x(0)=0  
    A_0(1)=1; 
    A_plus(1)=-1;
    A_minus(1)=0;
    d(1)=0; 
     

    %Second boundary condition v(L) = 0
    A_0(N)=1;
    A_plus(N)=0;
    A_minus(N)=0;
    d(N)=0;

    % future velocity prediction
    v_f_cnm=solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

%--------------------------------------------------------------------------
% Predict future v using explicit method
%--------------------------------------------------------------------------
function v_f = v_step_explicit(h_n,v_n,dx,dt,Oh)
    N = length(h_n);
    hx_n = D(h_n,dx);
    k_n = kappa(h_n,hx_n,dx,Oh);
    kx_n=D(k_n,dx);
    vx_n = D(v_n,dx);
    vxx_n = D2(v_n,dx);
   
    %future velcity
    A_n = (hx_n./h_n-v_n+v_n(1));
    dv_dt = A_n.*vx_n+vxx_n+kx_n;
    v_f = v_n+dt*dv_dt; 

    %boundary conditions 
    v_f(1)=v_f(2); %v_x_tip=0
    v_f(N)=0; % v_wall=0
end

%--------------------------------------------------------------------------
% Calculate curvature function (dimensionless definition)
%--------------------------------------------------------------------------
function k = kappa(h,hx,dx,Oh)
    hxx = D2(h,dx);
    x_0 = 4*Oh;
    k = hxx./(1+hx.^2/x_0^2).^(3/2)/x_0^2; % curvature
    
end

%--------------------------------------------------------------------------
% Predict Future h
%--------------------------------------------------------------------------
function h_f = h_step(h_n,v_n,v_f,dx,dt)
    N = length(h_n);
    A_0 = ones(N,1)'; %A(i,i) array 
    A_plus = zeros(N,1)';%A(i,i+1) array
    A_minus = zeros(N,1)';%A(i,i-1) array
    
    PDE_plus = dt/dx/4*(v_f(4:N)-v_f(1));
    PDE_minus = -dt/dx/4*(v_f(2:N-2)-v_f(1));
    
    A_plus(3:N-1)=PDE_plus;
    A_minus(3:N-1)= PDE_minus;
    %vector
    hx_n = D(h_n,dx);
    vx_n = D(v_n,dx);
    d = h_n - dt/2*(hx_n.*(v_n-v_n(1))+vx_n.*h_n); %explicit prediction for t+dt/2


    %enforce boundary conditions for v on the Matrix

    %first boundary condition h(0)=0
    A_0(1)=1;
    A_plus(1)=0;
    A_minus(1)=0;
    d(1)=0;
    
    %second boundary condition h(1)-h(2)/sqrt(2)=0 or h->Ax^(1/2) x->0
    
    A_0(2)=1;
    A_plus(2)=-1/sqrt(2);
    A_minus(2)=0;
    d(2)=0;
    

    %third boundary condition h(L) = 1
    
    A_0(N)=1;
    A_plus(N)=0;
    A_minus(N)=0;

    d(N)=h_n(N);
    
    h_f = solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

%--------------------------------------------------------------------------
% Solve Tridiagonal Matrix System
% cite: L.H. Thomas. Elliptic problems in linear difference equations
% over a network. Watson Sc. Comp. Lab. Rep., 1949.
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
% Return first partial derivative wrt position, central limit 
%--------------------------------------------------------------------------
function y=D(f,dx)
    N=length(f);
    df=linspace(0,0,N);
    for i = 2:N-1
        df(i)=(f(i+1)-f(i-1))/2;
    end
    y=df./dx;
    y(1)=2*y(2)-y(3);
    y(N)=2*y(N-1)-y(N-2);
end 

%--------------------------------------------------------------------------
% Return second partial derivative wrt position, central limit  
%--------------------------------------------------------------------------
function y=D2(f,dx)
    N=length(f);
    y=linspace(0,0,N);
    for i = 2:N-1
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
function plot_heights(t,x,h,Oh)  
    cmap = colororder(); %color sequence of time snaps
    Nframe = length(t);
    Ts = round([1,linspace(1,4,4)*Nframe/4]); % timesnaps showed in figure
    figure 
    for i = 1:length(Ts)
        label_name = 'T = '+string(round(t(Ts(i)),2));
        if i ==1
            label_name = sprintf('T = %2.1e',t(Ts(1)));
        end 
        plot(x,h(:,Ts(i)),'color', cmap(i,:), ...
            'linewidth',2, DisplayName=label_name);
        hold on
        plot(x,-h(:,Ts(i)),'color', cmap(i,:),'linewidth',2, ...
            'HandleVisibility','off');
        hold on
    end 
    label_figure('Position in Tip Frame, X''', 'Height, H(X'',T)',Oh);
    hold off
end


%--------------------------------------------------------------------------
% Plot Velocities
%--------------------------------------------------------------------------
function plot_velocities(t,x,v,Oh) 
    Nframe = length(t);
    Ts = round([1,linspace(1,4,4)*Nframe/4]); % timesnaps showed in figure 
    figure 
    for i = 1:length(Ts)
        label_name = 'T = '+string(round(t(Ts(i)),2));
        if i ==1
            label_name = sprintf('T = %2.1e',t(Ts(1)));
        end 
        
        plot(x,-v(:,Ts(i))+v(1,Ts(i)),'linewidth',2, ...
            DisplayName=label_name);
        hold on
    end     
    label_figure('Position in Tip Frame, X''', 'Velocity in Tip Frame, -V''(X'',T)',Oh);
    hold off
end

%--------------------------------------------------------------------------
% Plot Velocity of Free End Point
%--------------------------------------------------------------------------
function plot_tip_velocity(t,v,Oh)
    % neumerical plot
    figure
    plot(t,v(1,:),'linewidth', 2, DisplayName='Numerical V_{tip}'); 
    label_figure('time, T', 'Tip Velocity, U(T)',Oh);
    hold off
end


%--------------------------------------------------------------------------
% Plot Maximum Height
%--------------------------------------------------------------------------
function plot_maximum_height(t,h,Oh)     
    % neumerical plot
    hmax = linspace(0,0,length(t));
    for it=1:length(t)
        hmax(it)=max(h(:,it));
    end
    figure
    plot(t,hmax,'linewidth', 2, DisplayName='Numerical'); 
    
    label_figure('time, T ', 'Maximum Height, H_{max}(T)',Oh);
    hold off
end

%--------------------------------------------------------------------------
% Labeling the Figures
%--------------------------------------------------------------------------
function label_figure(x_label,y_label,Oh)
    % title of the calculations and figures
    fontsize = 18; %fontsize 
    %title_string = 'Oh = '+string(Oh);

    xlabel(x_label,'fontsize',fontsize,'fontname','times') %x-axis label
    ylabel(y_label,'fontsize',fontsize,'fontname','times') %y-axis label

    %title(title_string,'fontsize',fontsize,'Interpreter','latex') %titile
    legend('Location','best','Fontsize',fontsize,'fontname','times')

    grid on
    grid minor
    set(gca,'FontSize',fontsize,'fontname','times')
end


