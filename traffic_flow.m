x = linspace(0,1,1000);
t = linspace(0,100,10000);
f = simulate(x,t);
X = x(1:round(length(x)/100):end);
T = t(1:round(length(t)/100):end);
F = f(1:round(length(x)/100):end,1:round(length(t)/100):end)';

figure
[X,T]=meshgrid(X,T);
surf(X,T,F)
colormap winter
%shading interp
ax = gca;
ax.FontSize = 12; 
xlabel('x-axis','FontSize',15)
ylabel('time, t','FontSize',15)
zlabel('Function, f','FontSize',15)
grid on
grid minor

figure
for i = 1:round(length(t)/5):length(t)
    plot(x,f(:,i),'linewidth',2,'DisplayName','Time, t = '+string(round(t(i))))
    hold on
end
plot(x,f(:,length(t)),'linewidth',2,'DisplayName','Time, t = '+string(round(t(length(t)))))
hold off
legend('Location','best','Fontsize',18)
xlabel('x-axis','FontSize',18)
ylabel('Function, f','FontSize',18)
ax = gca;
ax.FontSize = 15; 
grid on
grid minor

figure
N = length(t);
cut = N/4;
plot(t,f(length(x),:), 'linewidth',2,'DisplayName','Numerical Simulation')
hold on
plot(t(cut:end),0.5./t(cut:end), 'linewidth',2, 'DisplayName','0.5/t')
hold off
xlabel('time, t','FontSize',18)
ylabel('Function, f(x=1,t)','FontSize',18)
legend('Location','best','Fontsize',18)

ax = gca;
ax.FontSize = 15; 
grid on
grid minor

function y=k(f)
    y = 2*f-3*f.^2;
end

function f = simulate(x,t)
    N = length(x);
    NT = length(t);
    dx = x(N)/N;
    dt = t(N)/N;
    f = zeros(N,NT);
    r = 50;
    f(:,1)=0.01*exp(-10*(x-1).^2);
    for it = 1:NT-1
        f(:,it+1)=step(f(:,it),N,dx,dt);
    end
end

function y = step_exp(f,dx,dt)
    y = f-dt/dx*k(f).*(D(f,dx));
end

function y = step(f,N,dx,dt)

    f_exp = step_exp(f,dx,dt);
    A_0 = zeros(N,1);
    A_plus = zeros(N,1);
    A_minus = zeros(N,1);
    A_0(1)=1+dt/dx/2*k(f_exp(1));
    A_plus(1)=-dt/dx/2*k(f_exp(1));
    for i = 2:N-1
        A_0(i)=1;
        A_plus(i)=k(f_exp(i))*dt/dx/4;
        A_minus(i)=-k(f_exp(i))*dt/dx/4;
    end
    A_0(N)=1+dt/dx/2*k(f_exp(N));
    A_minus(N)=-dt/dx/2*k(f_exp(N));

    d = f-1/2*dt/dx*k(f).*(D(f,dx));
    A_0(1)=1;
    A_minus(1)=-1;
    d(1)=0;
    
    y = solve_tridiagnoal(A_0,A_minus,A_plus,d);
end

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
% cite: L.H. Thomas. Elliptic problems in linear difference equations over a
% network. Watson Sc. Comp. Lab. Rep., 1949.
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
