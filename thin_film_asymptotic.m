N = 500; % number of space bins
L = 10; %initial length of film

dx = L/N; %bin gap
x = linspace(0,L,N);
NT = 1e6;%number of time bins
T = 100; %total time
dt = T/NT; %time gap

Nframe = min(NT,1e5); %number of frames to be saved

[t,h,u] = simulate(T,dt,dx,L,Nframe);

Ts = round(linspace(1,5,5)*length(t)/5); %more showed
Ts_vel = round(linspace(1,500,5)); %more showed
fsize = 18; %fontsize 

% title of the calculations and figures
s = 'L = %2.0f, Total Time = %2.0f, \\Delta x = %1.0e, \\Delta t = %1.0e';
title_string = sprintf(s,L,T,dx,dt);
plot_heights(t,x,h,Ts,fsize)
plot_velocities(t,x,h,Ts,fsize)   
plot_V_TC(t,u,fsize)  


function [t_save,h_save,u_save] = simulate(T,dt,dx,L,Nframe)
    N = round(L/dx);
    N_time = round(T/dt); 

    % the information will be saved here
    % we only save a fraction of the time data to save space
    % N_frame (No of saved time bins) < N_time (No of time bins)
    frame_step = round(N_time/Nframe); % time bins between two saved frame
    h_save = zeros(N,Nframe);
    t_save = linspace(0,T,Nframe);
    u_save = zeros(Nframe,1);

    h = linspace(0,0,N);
    
    t = 0;
    u = 0;
    for it = 0:N_time-1
        % save everything here
        if mod(it,frame_step) == 0
            h_save(:,round(it/frame_step)+1)=h;
            u_save(round(it/frame_step)+1)=u;
        end
        
        [hf,u] = h_cnm_2(h,t,dt,dx);

        % future data becomes present data
        h = hf;
        t = t+dt;
        
    end
end

function y = U(h,t,dx)
    %y = (3*t+h(3)-4*h(2))/(2*(1+t)*dx);
    hx = (h(2)-t)/dx;
    c = (dx-2*hx)/(1+t);
    y = (1-sqrt(1-c*dx))/dx;
end

function y = h_t(h,t,dx)
    y = U(h,t,dx)*D(h,dx)+D2(h,dx);
end

function y = h_exp(h,ht,t,dt)
    y = h+dt*ht;
    y(1)=t;
    y(length(y)) = 0;
end

function [y,d] = h_cnm_1(h,t,dt,dx)
    N = length(h);
    ht = h_t(h,t,dx);
    hf = h_exp(h,ht,t,dt);

    d = h_exp(h,ht,t,dt/2);
    d(1) = t+dt;
    d(N) = 0;

    u = U(hf,t+dt,dx);
    [A_0,A_plus,A_minus]= Matrix(dx,dt,u,N);
    y = solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

function [y,u] = h_cnm_2(h,t,dt,dx)
    N = length(h);
    [hf,d] = h_cnm_1(h,t,dt,dx);
    u = U(hf,t+dt,dx);
    [A_0,A_plus,A_minus]= Matrix(dx,dt,u,N);
    y = solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

function [A_0,A_plus,A_minus]= Matrix(dx,dt,u,N)

    A_0 = zeros(N,1)+1+dt/dx^2;
    A_plus = zeros(N,1)-dt/dx*(1/dx+u/2)/2;
    A_minus = zeros(N,1)-dt/dx*(1/dx-u/2)/2;

    A_0(1)=1;
    A_0(N)=1;
    
    A_plus(1)=0;
    A_plus(N)=0;

    A_minus(1)=0;
    A_minus(N)=0;
end

%--------------------------------------------------------------------------
% Solve Tridiagonal Matrix System
% cite: https://github.com/tamaskis/tridiagonal-MATLAB/blob/main/
% Tridiagonal_Matrix_Algorithm.pdf
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
    y=df/dx;
    y(1)=2*y(2)-y(3);
    y(N)=2*y(N-1)-y(N-2);
end 

%--------------------------------------------------------------------------
% Return second partial derivative wrt position, central limit  
%--------------------------------------------------------------------------
function y=D2(f,dx)
    N=length(f);
    d2f=linspace(0,0,N);
    for i = 2:N-1
        d2f(i)=f(i+1)+f(i-1)-2*f(i);  
    end
    y=d2f/dx^2;
    y(1)=2*y(2)-y(3);
    y(N)=2*y(N-1)-y(N-2);
end 

function plot_heights(t,x,h,Ts,fsize)  
    figure 
    for i = 1:length(Ts)
        plot(x,1+h(:,Ts(i)), 'linewidth',1.5,DisplayName='T = '+string(round(t(Ts(i)),2)));
        hold on
    end
    xp = x(1:10:end);
    T = t(length(t));
    plot(xp,1+T*exp(-xp),'o','linewidth',1.5,DisplayName = '1+Te^{-X''}')
    hold off
    xlim([0,10])
    ylim([1,t(length(t))+1])
    xlabel('Position in Moving Tip Frame, X''','fontsize',fsize,'fontname','times')
    ylabel('Height, H','fontsize',fsize,'fontname','times')
    legend('location','best','fontsize',fsize,'fontname','times')
    grid on
    grid minor
    set(gca,'FontSize',fsize,'fontname','times')
end

function plot_velocities(t,x,h,Ts,fsize)  
    dx = x(2)-x(1);
    figure 
    for i = 1:length(Ts)
        H = h(:,Ts(i));
        Ti = t(Ts(i));
        u = -U(H,Ti,dx);
        v = u-D(H,dx)./(H+1)';
        v(length(v))=u;
        plot(x,-v, 'linewidth',1.5,DisplayName='T = '+string(round(Ti,2)));
        hold on
    end
    xp = x(1:10:end);
    T = t(length(t));
    plot(xp,1./(T*exp(-xp)+1),'o','linewidth',1.5,DisplayName = '1/(Te^{-X''}+1)')
    hold off
    xlim([0,10])
    xlabel('Position in Moving Tip Frame, X''','fontsize',fsize,'fontname','times')
    ylabel('Velocity in Moving Tip Frame, -V''','fontsize',fsize,'fontname','times')
    legend('location','best','fontsize',fsize,'fontname','times')
    grid on
    grid minor
    set(gca,'FontSize',fsize,'fontname','times')
end

function plot_V_TC(t,u,fsize)  
    figure 
    loglog(t,u,'color', 'black', 'linewidth',1.5)
    hold on
    tp = 10.^(linspace(-3,0,20));
    loglog(tp,(4*tp/pi).^(1/2),'o','linewidth',2)
    hold off
    ylabel('Tip Velocity, U(T)','fontsize',fsize,'fontname','times')
    xlabel('time, T','fontsize',fsize,'fontname','times')
    legend('Numerical', '\surd(4T/\pi)','location','southeast','fontsize',fsize,'fontname','times')
    grid on
    grid minor
    set(gca,'FontSize',fsize,'fontname','times')
end
