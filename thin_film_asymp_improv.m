L = 10;
T = 15;
N = 500;
NT = 1e5;

run_code(L,T,N,NT)

% Run this function to simulate and plot 
% L is our l_0/h_0/4/Oh
% T is total time
% N is number of spatial bins
% NT is number of temporal bins
function run_code(L,T,N,NT)

    Nframe = min(NT,10^5);
    tic
    [t,x,h,v]=simulate(L,N,T,NT,Nframe);
    toc
    plot_graph(x,h,T,Nframe,L,2,1)
    plot_graph(x,v,T,Nframe,L,2,2)
    plot_veltip(t,v,N,L,2)
    %plot_hwall(t,h,L,2
    %plot_length(t,x,N,L,2)
end

% Run this function to simulate and plot
% for multiple values of L, T, N, NT. set them in an array
% plot all = 1 will plot the graphs for all of them
% plot all != 1 will only plot the maximum velocity reached vs L 
function run_code_array(L,T,N,NT,plot_all)
    transition = zeros(length(L),1);
    v_max = zeros(length(L),1);
    t_max = zeros(length(L),1);    
    for i = 1:length(L)
        Nframe = min(NT(i),10^5);
        [t,x,h,v]=simulate(L(i),N(i),T(i),NT(i),Nframe);
        if plot_all == 1
           plot_graph(x,h,T(i),Nframe,L(i),2,'Height')
           plot_graph(x,v,T(i),Nframe,L(i),2,'Velocity')
           plot_hwall(t,h,L(i),2)
           plot_veltip(t,v,N(i),L(i),2)
           plot_length(t,x,N(i),L(i),2)
        end
        [v_max(i),imax]=max(v(N(i),:));
        t_max(i)=t(imax);
        dT = t(2)-t(1); 
        phase=D(v(N,:)',dT);
        [~,itrans]=min(phase);
        transition(i) = t(itrans);
    
    end
    plot_transition(L,transition)
    plot_vmax(L,v_max)
    plot_tmax(L,t_max)    
end

function [t,x,h,v]=simulate(L,N,T,NT,Nframe)
    h = zeros(N,Nframe);
    x = zeros(N,Nframe);
    v = zeros(N,Nframe);
    t = linspace(0,T,Nframe);

    frame_step = round(NT/Nframe);
    
    sn = L^2;
    hn = linspace(1,1,N);
    x0 = linspace(0,1,N);
    tn = 0;

    dt = T/NT;
    
    iframe = 1;

    for it = 1:NT   
        if mod(it,frame_step) == 0
            h(:,iframe)=hn*(1+tn+dt);
            x(:,iframe)=x0*sqrt(sn);
            v(:,iframe)=v_calculate(h(:,iframe),x(:,iframe));
            iframe = iframe+1;
        end
        [hn,sn]=step(dt,x0,sn,tn,hn);
        tn = tn + dt;
    end  
end

function [hf,sf] = step(dt,x,sn,tn,hn)
    dx = x(2)-x(1);

    sf = s_step(dt,dx,sn,tn,hn,hn); % rough update of sf with hf = hn
    hf = h_step(dt,x,sn,sf,tn,hn); % rough update hf
    
    sf = s_step(dt,dx,sn,tn,hn,hf); % accurate update of sf
    hf = h_step(dt,x,sn,sf,tn,hn); % accurate update hf
end

function sf = s_step(dt,dx,sn,tn,hn,hf)
    N = length(hn);
    hx = (1-hn(N-1)/2-hf(N-1)/2)/dx;
    t = tn + dt/2;

    c = 2*hx + sn*dx/(1+t);
    b = 1+dx*dt/(1+t)/2;
    a = dx/4;

    st = (-b+sqrt(b^2-4*a*c))/(2*a); 
    sf = sn + st*dt;
end

function hf = h_step(dt,x,sn,sf,tn,hn)
    dx = x(2)-x(1);
    N = length(x);

    t = tn+dt/2;
    st = (sf-sn)/dt;
    s = (sf+sn)/2;
    
    R = dt/dx^2/s;
    H = dt/(1+t)/2;
    G = st*dt*x/(8*s*dx);
    
    A_0 = [1,(1+R+H)*ones(1,N-2),1];
    A_plus = [-1,-(R/2+G(2:N-1)),0];
    A_minus = [0,-(R/2-G(2:N-1)),0];

    ht = 1/s*D2(hn,dx)-hn/(1+t)+st*x.*D(hn,dx)/(2*s);

    d = [0,hn(2:N-1) + dt/2*ht(2:N-1),1];

    hf = solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

function v = v_calculate(h,x)
    dx = x(2)-x(1);
    v = (D(log(h),dx))';
end

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

function y = D(f,dx)
    N=length(f);
    df=linspace(0,0,N);
    for i = 2:N-1
        df(i)=(f(i+1)-f(i-1))/2;
    end
    y=df/dx;
    y3_2 = (f(2)-f(1))/dx;
    yN = (f(N)-f(N-1))/dx;
    y(1)=2*y3_2-y(2);
    y(N)=2*yN-y(N-1);
end

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

function plot_graph(x,h,T,Nframe,L,lw,ylab)
    if ylab == 1
        label = 'Height';
    else
        label = 'Velocity';
    end
    N = length(x(:,1));
    figure 
    plot(x(:,1),h(:,1),DisplayName = 'time = '+string(T/Nframe),LineWidth=lw)
    hold on
    for i = round(Nframe*0.2):round(Nframe*0.2):Nframe
        plot(x(:,i),h(:,i),DisplayName = 'time = '+string(i*T/Nframe),LineWidth=lw)
        hold on
    end

    if (L >= 10) && (T<L)   
        Lt = x(N,Nframe);
        xplot = linspace(0,Lt,50);
        xp = Lt-xplot;

        if ylab == 1
            hcomp = T*exp(-(T+1)*xp/T)+1;
            plot(xplot,hcomp,'o',DisplayName = 'Analytic, time = '+string(T));
        else
            vcomp = (T+1)*exp(-(T+1)*xp/T)./(T*exp(-(T+1)*xp/T)+1);
            plot(xplot,vcomp,'o',DisplayName = 'Analytic, time = '+string(T));
        end
        hold off
    end

    ax = gca;
    ax.FontSize = 15; 
    l=legend(FontSize=15,Location="best");
    set(l,'Interpreter','Latex');
    xlabel('x-axis',FontSize=15)
    ylabel(label,FontSize=15) 
    title('$\mathcal{L} = $ '+string(L),'Interpreter','latex',FontSize=15)
    grid on
    grid minor
end

function plot_veltip(t,v,N,L,lw)
    NT=length(t);
    figure
    loglog(t,v(N,:),DisplayName = 'Numerical',LineWidth=lw)
    hold on
    loglog(t,sqrt(4*t/pi), '--',DisplayName='$\sqrt{\frac{4T}{\pi}}$',LineWidth=lw)
    hold on
    loglog(t,L./(1+t).^2, '-.',DisplayName='$\frac{\mathcal{L}}{(1+T)^2}$',LineWidth=lw)
    hold on
    loglog(t,1-1./(t+1).^2, ':',DisplayName='$1-\frac{1}{(1+T)^2}$',LineWidth=lw)
    hold off
    vmin = min(v(N,2:end));
    vmax = max(v(N,2:end));
    ylim([vmin,vmax])
    xlim([t(2),t(NT)])
    ax = gca;
    ax.FontSize = 15; 
    l=legend(FontSize=15,Location="best");
    set(l,'Interpreter','Latex');
    xlabel('time, T',FontSize=15)
    ylabel('Tip Velocity, U(T)',FontSize=15)
    title('$\mathcal{L} = $ '+string(L),'Interpreter','latex',FontSize=15)
    grid on
    grid minor
end

function plot_length(t,x,N,L,lw)
    figure
    plot(t,L-x(N,:),LineWidth=lw)
    hold off
    
    ax = gca;
    ax.FontSize = 15; 
    xlabel('time, T',FontSize=15)
    ylabel('$\Delta \mathcal{L}$','Interpreter','latex',FontSize=15)
    title('$\mathcal{L} = $ '+string(L),'Interpreter','latex',FontSize=15)
    grid on
    grid minor
end

function plot_hwall(t,h,L,lw)
    figure
    plot(t,h(1,:),LineWidth=lw,DisplayName='Numerical')
    hold on
    plot(t,1+t,'--',LineWidth=lw,DisplayName='1+T')
    hold off

    l=legend(FontSize=15,Location="best");
    set(l,'Interpreter','Latex');
    ax = gca;
    ax.FontSize = 15; 
    xlabel('time, T',FontSize=15)
    ylabel('Height at Wall',FontSize=15)
    title('$\mathcal{L} = $ '+string(L),'Interpreter','latex',FontSize=15)
    grid on
    grid minor
end

function plot_vmax(L,v_max)
    figure
    loglog(L,v_max,'o',LineWidth=2,DisplayName='Simulation')
    hold on
    loglog(L,L/4,'--',LineWidth=2,DisplayName='Analytic, $\mathcal{L}\ll1$')
    hold on
    loglog(L,L.^0,':',LineWidth=2, DisplayName='Taylor-Culick Velocity')
    hold off
    ylim([0,1.1])
    l=legend(FontSize=15,Location="best");
    set(l,'Interpreter','Latex');
    ax = gca;
    ax.FontSize = 15; 
    xlabel('$\mathcal{L}$ ','Interpreter','latex',FontSize=15)
    ylabel('Maximum Velocity',FontSize=15)
    grid on
    grid minor
end

function plot_tmax(L,t_max)  
    figure
    loglog(L,t_max,'o',LineWidth=2,DisplayName='Simulation')
    hold on
    loglog(L,L,'--',LineWidth=2, DisplayName='Analytic, $\mathcal{L}\gg1$')
    hold off
    ylim([0,max(t_max)])
    l=legend(FontSize=15,Location="best");
    set(l,'Interpreter','Latex');
    ax = gca;
    ax.FontSize = 15; 
    xlabel('$\mathcal{L}$ ','Interpreter','latex',FontSize=15)
    ylabel('Time to reach Maximum Velocity',FontSize=15)
    grid on
    grid minor
end

function plot_transition(L,transition)
    figure
    loglog(L,transition,'o',LineWidth=2,DisplayName='Simulation')
    hold on
    loglog(L,L,'--',LineWidth=2, DisplayName='Analytic, $\mathcal{L}\gg1$')
    hold off
    l=legend(FontSize=15,Location="best");
    set(l,'Interpreter','Latex');
    ax = gca;
    ax.FontSize = 15; 
    xlabel('$\mathcal{L}$ ','Interpreter','latex',FontSize=15)
    ylabel('Time of Phase Change',FontSize=15)
    grid on
    grid minor
end
