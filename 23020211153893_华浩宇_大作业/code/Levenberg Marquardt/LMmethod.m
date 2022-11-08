clear all;clc;close all;
syms c1 c2 c3 c4 t y
data=importdata('efit1.dat');
t1=data.data(:,1);
y1=data.data(:,2);
epsilon=1e-8;
n=200;

disp('Example 3.7')
c0=[-1;-2;1;-1];
lambda=0.001;
c_root1=LevenbergMarquardt(t1,y1,c0,epsilon,epsilon,lambda,n)

function [ x ] = LevenbergMarquardt( t , y , x0 , epsilon_1 , epsilon_2 , tau , N )
%% Levenberg-Marquardt算法，进行高斯曲线拟合
% 输入变量为 t , y , x0 , epsilon_1 , epsilon_2 , tau , N
% t：自变量
% y：因变量
% x0：参数初始猜测值
% epsilon_1：迭代终止误差
% epsilon_2：迭代终止误差
% tau：倍率
% N：迭代最高次数
% 输出变量为 x
% x：迭代终止后参数的结果
Fx = [];
gg = [];
mul = [];
%% Objective function vector
fun = @(x) y-(x(3)*exp(x(1)*t)+x(4)*exp(x(2)*t)); % gauss type
% fun = @(x) log( abs( y - x(4) ) ) - log( abs( x(1) ) ) + ( ( t - x(2) ) / x(3) ).^2 / 2; % log type
%% Jacobi matrix for gauss function
Jacobi = @(x) [
    -t*x(3).*exp(x(1)*t),...
    -t*x(4).*exp(x(2)*t),...
    -exp(x(1)*t),...
    -exp(x(2)*t)
    ]; % gauss type
% Jacobi = @(x) -1 * [
%     ones( length(t) , 1 ) / x(1) ,...
%     ( x(2) - t ) / x(3)^2,...
%     -( t - x(2) ).^2 / x(3)^3,...
%     1 ./ ( x(4) - y )
%     ]; % log type
%%
count = 0;
v = 2;
x = x0;
Jmatrix = Jacobi( x );
A = Jmatrix' * Jmatrix;
g = Jmatrix' * fun( x );
found = ( norm( g , Inf ) < epsilon_1 );
mu = tau * max( diag( A ) );

while ~found && ( count < N )
    count = count + 1;
    h_lm = ( A + mu * diag( ones( 4 , 1 ) ) ) \ ( -g );
    if ( norm( h_lm ) <= ( epsilon_2 * ( norm( x ) + epsilon_2 ) ) )
        found = 1;
    else
        x_new = x + h_lm;
        gain_radio = ( fun( x )' * fun( x ) - fun( x_new )' * fun( x_new ) ) / 2 /...
            ( 0.5 * h_lm' * ( mu * h_lm - g ) );
        if gain_radio > 0
            x = x_new;
            Jmatrix = Jacobi( x );
            A = Jmatrix' * Jmatrix;
            g = Jmatrix' * fun( x );
            found = ( norm( g , Inf ) <= epsilon_1 );
            gg = [gg norm(g,Inf)];
            mul = [mul mu];
            tmp = 0.5*fun(x)'*fun(x);
            Fx = [Fx tmp];
            mu = mu * max( [ 1/3 , 1 - ( 2 * gain_radio - 1 )^3 ] );
            v = 2;

        else
            mu = mu * v;
            v = 2 * v;
        end
    end
end

figure(1);
hold on;
plot(t,y,'r.');
ax = 0:0.01:1;
plot(ax,x(3)*exp(x(1)*ax)+x(4)*exp(x(2)*ax),'b');
plot(ax,x0(3)*exp(x0(1)*ax)+x0(4)*exp(x0(2)*ax),'k');
hold off;

figure(2);
hold on;
plot(1:1:count,Fx,'bp')
plot(1:1:count,gg,'rv')
plot(1:1:count,mul,'ko')
hold off;
end
