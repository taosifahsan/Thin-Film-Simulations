% create matrix A for implicit part
% these three digonal arrays are only non zero entries of Matrix A
%{
repeat = 100;
cut = 6;
start_value = 18;
end_value = 41;
sample_no = (end_value-start_value)+1;
m = linspace(start_value,end_value,sample_no);


N_array = round(1.25.^m);
runtime_direct = linspace(0,0,length(m));
runtime_error_direct = linspace(0,0,length(m));
runtime_tri = linspace(0,0,length(m));
runtime_error_tri = linspace(0,0,length(m));
for i = 1:length(m)
    N = N_array(i);
    range = 10^10;
    
    A = zeros(N,N);

    A_0 = linspace(0,0,N); %A(i,i) array
    A_plus = linspace(0,0,N); %A(i,i+1) array
    A_minus = linspace(0,0,N); %A(i,i-1) array


    time_direct = linspace(0,0,repeat);
    time_tri = linspace(0,0,repeat);

    for j = 1:repeat
        d = 2*range*(rand(N,1))'-1;
        for ix = 2:N-1
            A(ix,ix)=2*range*rand(1)-1;
            A(ix,ix+1)=2*range*rand(1)-1;
            A(ix,ix-1)=2*range*rand(1)-1;
    
            A_0(ix)=A(ix,ix); 
            A_plus(ix)=A(ix,ix+1);  
            A_minus(ix)=A(ix,ix-1);
        end   
    
        A(1,1)=2*range*rand(1)-1;
        A(N,N-1)=2*range*rand(1)-1;
        A(N,N)=2*range*rand(1)-1;
    
        A_0(1)=A(1,1); 
        A_minus(N)=A(N,N-1);
        A_0(N)=A(N,N);

        tic 
        x1 = (A\d')';
        time_direct(j) = toc;
        
        tic 
        x2 = solve_tridiagnoal(A_0,A_minus,A_plus,d);
        time_tri(j) = toc;
    end

    runtime_direct(i) = mean(time_direct);
    runtime_tri(i) = mean(time_tri);

    runtime_error_direct(i) = sqrt(var(time_direct));
    runtime_error_tri(i) = sqrt(var(time_tri));
end

Y_direct = log(runtime_direct)/log(10);
err_direct = log(1+runtime_error_direct./runtime_direct)/log(10);
Y_tri = log(runtime_tri)/log(10);
err_tri = log(1+runtime_error_tri./runtime_tri)/log(10);
X = log(N_array)/log(10);

coeff_direct=polyfit(X(cut:end),Y_direct(cut:end),1);
str_direct = 'slope = '+string(round(coeff_direct(1),2));
coeff_tri=polyfit(X(cut:end),Y_tri(cut:end),1);
str_tri = 'slope = '+string(round(coeff_tri(1),2));

fit_direct = coeff_direct(2)+coeff_direct(1)*X;
fit_tri = coeff_tri(2)+coeff_tri(1)*X;
%}
fsize = 18;
figure
fig=errorbar(X,Y_direct, err_direct,'o','linewidth',1.5);
hold on
errorbar(X,Y_tri, err_tri,'o','linewidth',1.5);
hold on
plot(X,fit_direct,'linewidth',1.5,'Color','blue')
hold on
plot(X,fit_tri,'linewidth',1.5,'Color','red')
hold on 
plot([X(cut),X(cut)],[-8,1])
hold off
xlabel('log_{10}N', 'fontsize',fsize,'fontname','times')
ylabel('log_{10}(Run Time/s)', 'fontsize',fsize,'fontname','times')
legend('Direct Inversion', 'Tridiagonal Inversion', ...
    str_direct,str_tri,'Location','southeast','fontname','times')
set(gca, 'fontsize',fsize,'fontname','times')
title({'Number of Trials per point = '+string(repeat), ...
    'Cut off for fit, N > '+string(N_array(cut))}, 'fontsize',fsize,'fontname','times')
%xlim([0,N_array(length(N_array))])
grid on
grid minor
%saveas(fig,"Direct_vs_Tridiagonal_inversion_runtime.png")
%saveas(fig,"Direct_vs_Tridiagonal_inversion_runtime.fig",'fig')

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