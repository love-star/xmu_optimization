clear all
close all
clc

M=100;

data=importdata('efit1.dat');
t_data=data.data(:,1)';
y_data=data.data(:,2)';

x0=[-1 -2 1 -1];
k_max=200;
epsilon_1=1.0e-8;
epsilon_2=1.0e-8;
epsilon_3=1.0e-8;
delta=1;

disp('Example 3.11');
fit_result=powells_dog_leg_sinus(t_data, y_data, x0, k_max, epsilon_1, epsilon_2, epsilon_3, delta)

function [x_new] = powells_dog_leg_sinus(t_data, y_data, x0, k_max, epsilon_1, epsilon_2, epsilon_3, delta)

    M=length(t_data(1,:));
    k=0;
    x=x0;
    F_x_list = [];
    g_norm_list = [];
    delta_list = [];
    % Jacobian matrix
    J=zeros(M,4);
    J(:,1)= -t_data*x(3).*exp(x(1)*t_data);
    J(:,2)= -t_data*x(4).*exp(x(2)*t_data);
    J(:,3)= -exp(x(1)*t_data);
    J(:,4)= -exp(x(2)*t_data);
    % calculate f
    f=zeros(1,M);
    f=y_data-(x(3)*exp(x(1)*t_data)+x(4)*exp(x(2)*t_data));
    % calculate g
    g=transpose(J(:,:))*transpose(f(1,:));
    % calculate norm f
    f_norm=sqrt(sum(f(1,:).*f(1,:)));
    % calculate norm g
    g_norm=sqrt(sum(g(:,1).*g(:,1)));
    % calculate boolean variable found
    found_bool=(f_norm <= epsilon_3) | (g_norm <= epsilon_1);
    % main loop
    while ((~found_bool) && (k < k_max))
        % increase iteration variable k
        k=k+1;
        % Jacobian matrix
        J=zeros(M,4);
        % caclulate derivation according to amplitude
        % fill first column of Jacobian matrix
        J(:,1)= -t_data*x(3).*exp(x(1)*t_data);
        J(:,2)= -t_data*x(4).*exp(x(2)*t_data);
        J(:,3)= -exp(x(1)*t_data);
        J(:,4)= -exp(x(2)*t_data);
        % calculate f
        f=zeros(1,M);
        f=y_data-(x(3)*exp(x(1)*t_data)+x(4)*exp(x(2)*t_data));
        % calculate g
        g=transpose(J(:,:))*transpose(f(1,:));
        % calculate norm g
        g_norm=sqrt(sum(g(:,1).*g(:,1)));
        % calculate J(x)*g(x)
        Jg=J(:,:)*g(:,1);
        % calculate norm Jg
        Jg_norm=sqrt(sum(Jg(:,1).*Jg(:,1)));
        % calculate alpha
        alpha=g_norm^2/Jg_norm^2;
        % calculate h_sd
        h_sd=(-1.0)*alpha*g(:,1);
        h_sd=transpose(h_sd);
        % calculate h_gn
        A=transpose(J(:,:))*J(:,:);
        B=transpose(J(:,:))*transpose(f(1,:));
        h_gn=(-1.0)*transpose(B)*inv(A);
        % calculate h_dl
        % calculate norm h_gn
        h_gn_norm=sqrt(sum(h_gn(1,:).*h_gn(1,:)));
        % calculate norm alpha*h_sd
        alpha_h_sd=alpha*h_sd(1,:);
        alpha_h_sd_norm=sqrt(sum(alpha_h_sd(1,:).*alpha_h_sd(1,:)));
        % calculate norm h_sd
        h_sd_norm=sqrt(sum(h_sd(1,:).*h_sd(1,:)));
        % caclulate h_dl
        if (h_gn_norm <= delta)
            h_dl=h_gn(1,:);
        elseif (alpha_h_sd_norm >= delta)
            h_dl=(delta/h_sd_norm)*h_sd(1,:);
        else
            beta=sym('beta');
            h_dl=alpha*h_sd(1,:)+beta*(h_gn(1,:)-alpha*h_sd(1,:));
            h_dl_norm=sqrt(sum(h_dl(1,:).*h_dl(1,:)));
            eqn=(h_dl_norm == delta);
            beta_result=solve(eqn, beta, 'Real', true);
            beta_result=double(beta_result);
            clear beta;
            h_dl=alpha*h_sd(1,:)+beta_result(2,1)*(h_gn(1,:)-alpha*h_sd(1,:));
        end
        % calculate norm h_dl
        h_dl_norm=sqrt(sum(h_dl(1,:).*h_dl(1,:)));
        % calculate norm x
        x_norm=sqrt(sum(x(1,:).*x(1,:)));
        delta_list = [delta_list delta];
        g_norm_list = [g_norm_list g_norm];
        if (h_dl_norm <= epsilon_2*(x_norm+epsilon_2))
            found_bool=1;
        else
            x_new=x(1,:)+h_dl(1,:);
            % calculate F(x)
            % caclulate function f
            f=zeros(1,M);
            f=y_data-(x(3)*exp(x(1)*t_data)+x(4)*exp(x(2)*t_data));
            F_x=0.5*sum(f(1,:).*f(1,:));
            F_x_list = [F_x_list F_x];
            % calculate F(x_new)
            f=zeros(1,M);
            f=y_data-(x_new(3)*exp(x_new(1)*t_data)+x_new(4)*exp(x_new(2)*t_data));
            F_x_new=0.5*sum(f(1,:).*f(1,:));
            % ro denominator
            ro_denominator=0.5*h_dl(1,:)*(alpha*transpose(h_dl(1,:))-g(:,1));
            % calculate ro - gain ratio
            ro=(F_x-F_x_new)/ro_denominator;
            if (ro > 0)
                x=x_new(1,:);
                % Jacobian matrix
                J=zeros(M,4);
                % caclulate derivation according to amplitude
                % fill first column of Jacobian matrix
                J(:,1)= -t_data*x(3).*exp(x(1)*t_data);
                J(:,2)= -t_data*x(4).*exp(x(2)*t_data);
                J(:,3)= -exp(x(1)*t_data);
                J(:,4)= -exp(x(2)*t_data);
                % calculate f
                f=zeros(1,M);
                f=y_data-(x(3)*exp(x(1)*t_data)+x(4)*exp(x(2)*t_data));
                % calculate norm f
                f_norm=sqrt(sum(f(1,:).*f(1,:)));
                % calculate g
                g=transpose(J(:,:))*transpose(f(1,:));
                % calculate norm g
                g_norm=sqrt(sum(g(:,1).*g(:,1)));
                % cacluate boolean variable found
                found_bool=(f_norm <= epsilon_3) | (g_norm <= epsilon_1);
            end
            if (ro > 0.75)
                % calculate norm h_dl
                h_dl_norm=sqrt(sum(h_dl(1,:).*h_dl(1,:)));
                % calculate new delta
                delta=max([delta 3*h_dl_norm]);
            elseif (ro < 0.25)
                % calculate new delta
                delta=delta/2;
                % calculate norm x
                x_norm=sqrt(sum(x(1,:).*x(1,:)));
                % calculate boolean variable found
                found_bool=(delta <= epsilon_2*(x_norm+epsilon_2));
            end
        end
    end

figure(1);
hold on;
plot(t_data,y_data,'r.');
ax = 0:0.01:1;
plot(ax,x0(3)*exp(x0(1)*ax)+x0(4)*exp(x0(2)*ax),'b');
plot(ax,x_new(3)*exp(x_new(1)*ax)+x_new(4)*exp(x_new(2)*ax),'k');
hold off;

figure(2);
hold on;
plot(1:1:k,F_x_list,'bp')
plot(1:1:k,g_norm_list,'rv')
plot(1:1:k,delta_list,'ko')
hold off;
end
