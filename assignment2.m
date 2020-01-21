%%
%Q1, a
%Uniform U[-2pi,pi]
y = unia(1000000);
histogram(y, 100)
hold on
y = 10000;
line([-2*pi,pi],[y,y],'LineWidth',2)
title('Uniform distribution over the interval [?2\pi, \pi]')
xlabel('-2\pi < x < \pi') 
ylabel('Frequency')
legend('Sample','Expected')

%%
%Q1, b
y = unib(100000)
histogram(y,100)
hold on
y1 = 5000;
line([1,2],[y1,y1],'LineWidth',2)
line([3,4],[y1,y1],'LineWidth',2)
line([5,6],[y1,y1],'LineWidth',2)
title('Uniform distribution over the union of the three intervals [1, 2]?[3, 4]?[5, 6]')
xlabel('x')
ylabel('Frequency')
legend('Sample','Expected')

%%
%Q1, c
a = box_muller(12.5,3)
histfit(a,100,'normal')
title('Gaussian distribution')
xlabel('x')
ylabel('Frequency')
legend('Sample','Expected')

%%
%Q1, d
a = genpoisson(0.7)
histfit(a,100,'Exponential')
title('Exponential distribution, lambda = 0.7')
xlabel('x')
ylabel('Frequency')
legend('Sample','Expected')

%%
a = genpoisson(1.5)
histfit(a,100,'Exponential')
title('Exponential distribution, lambda = 1.5')
xlabel('x')
ylabel('Frequency')
legend('Sample','Expected')

%%
a = genpoisson(3.5)
histfit(a,100,'Exponential')
title('Exponential distribution, lambda = 3.5')
xlabel('x')
ylabel('Frequency')
legend('Sample','Expected')

%%
%Q1, e
clear all;
for i = 1:10000
    a(i) = -1/2 + sqrt(6*rand+1/4);
end
histogram(a, 100)
hold on
x=0:2/10000:2
y=100/3+200/3*x
plot(x,y,'LineWidth',2)
title('Distribution over the interval [0,2]')
xlabel('0 ? x ? 2')
ylabel('Frequency')
legend('Sample','Expected')

%%
%Q2, a
%Check out function 'accreja' at the bottom of this file.

%%
%Q2, b
a = accreja(10000)
yyaxis left
histogram(a, 100)
ylabel('Frequency')
hold on
yyaxis right
fplot(@(x) (100+4/5*x^3-6*x^2)/1000, [0,10], 'LineWidth', 2)
ylim([0 0.3])
title('Distribution over the interval [0,10]')
xlabel('0 ? x ? 10') 
ylabel('f(x)')
legend('Sample','Expected')

%%
%Q2, c
%Bottom of this file is the function (accrejc).

%%
%Q2, d
[a,rate] = accrejc(10000)
yyaxis left
histogram(a, 100)
ylabel('Frequency')
hold on
yyaxis right
fplot(@(x) (100+4/5*x^3-6*x^2)/1000, [0,10], 'LineWidth', 2)
ylim([0 0.35])
title('Distribution over the interval [0,10]')
xlabel('0 ? x ? 10') 
ylabel('f(x)')
legend('Sample','Expected')
disp(rate)

%%
%Q3, a
[ans,se] = meanvalue(1000)

%%
%Q3, b
[ans,se] = importance_sampling(1000)

%%
%Q3, c
[ans,se] = isf1(1000)

%%
[ans,se] = isftwo(1000)

%%
[ans,se] = isfthree(1000)

%%
%Q4, a
%(i)
y = [0:1:100];
hold on
for i = 1:10
    x = rw(0.5,101);
    plot(y,x)
end
title('10 random walks, p = 0.5')
xlabel('Time step i')
ylabel('X(i)')

%%
%(ii)
y = zeros(100,1)
for i = 1:100
    x = rw(0.5,6)
    y(i) = x(6)
end
edges = [-5.5:1:5.5];
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-5:2:5];
j = 1
sample = zeros(6,1)
for i = -5:2:5
    sample(j) = bindist(i,5) * 100;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 100 5-step-walk end up')
xlabel('X(5)')
legend('Sampled','Expected')

%%
%(iii)
y = zeros(100,1)
for i = 1:100
    x = rw(0.5,21)
    y(i) = x(21)
end
edges = [-20.5:1:20.5]
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-20:2:20];
j = 1;
sample = zeros(21,1)
for i = -20:2:20
    sample(j) = bindist(i,20) * 100;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 100 20-step-walk end up')
xlabel('X_{20}')
legend('Sampled','Expected')

%%
%(iv)
y = zeros(1000,1)
for i = 1:1000
    x = rw(0.5,6)
    y(i) = x(6)
end
edges = [-5.5:1:5.5];
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-5:2:5];
j = 1
sample = zeros(6,1)
for i = -5:2:5
    sample(j) = bindist(i,5) * 1000;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 1000 5-step-walk end up')
xlabel('X(5)')
legend('Sampled','Expected')

%%
y = zeros(1000,1)
for i = 1:1000
    x = rw(0.5,21)
    y(i) = x(21)
end
edges = [-20.5:1:20.5]
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-20:2:20];
j = 1;
sample = zeros(21,1)
for i = -20:2:20
    sample(j) = bindist(i,20) * 1000;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 1000 20-step-walk end up')
xlabel('X_{20}')
legend('Sampled','Expected')

%%
%Q4, b
%(i)
y = [0:1:100];
hold on
for i = 1:100
    x = rw(0.2,101);
    plot(y,x)
end
title('100 random walks, p = 0.2')
xlabel('Time step i')
ylabel('X(i)')

%%
%(ii)
%p = 0.2
y = zeros(100,1)
for i = 1:100
    x = rw(0.2,6)
    y(i) = x(6)
end
edges = [-5.5:1:5.5]
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-5:2:5];
j = 1;
sample = zeros(6,1)
for i = -5:2:5
    sample(j) = pbindist(i,5,0.2) * 100;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 100 5-step-walk end up, p = 0.2')
xlabel('X_{5}')
legend('Sampled','Expected')

%%
%p = 0.1
y = zeros(100,1)
for i = 1:100
    x = rw(0.1,6)
    y(i) = x(6)
end
edges = [-5.5:1:5.5]
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-5:2:5];
j = 1;
sample = zeros(6,1)
for i = -5:2:5
    sample(j) = pbindist(i,5,0.1) * 100;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 100 5-step-walk end up, p = 0.1')
xlabel('X_{5}')
legend('Sampled','Expected')

%%
%p = 0.1
y = zeros(100,1)
for i = 1:100
    x = rw(0.8,6)
    y(i) = x(6)
end
edges = [-5.5:1:5.5]
histogram(y, edges)
ylabel('Frequency')
hold on
%Generate an expected sample:
k = [-5:2:5];
j = 1;
sample = zeros(6,1)
for i = -5:2:5
    sample(j) = pbindist(i,5,0.8) * 100;
    j = j + 1;
end
plot(k, sample, 'r', 'linewidth', 4);
title('Where 100 5-step-walk end up, p = 0.8')
xlabel('X_{5}')
legend('Sampled','Expected')

%%
%Functions:
function x = unia(N)
    x = zeros(1,N);
    for i = 1:N
        x(i) = 3*pi*rand - 2*pi;
    end
end

function x = unib(N)
    N = N*3
    x = zeros(1,N);
    for i = 1:N/3
        x(i) = rand + 1;
    end
    for i = N/3+1:N/3*2
        x(i) = rand + 3;
    end
    for i = N/3*2+1:N
        x(i) = rand + 5;
    end
end

function a = box_muller(mu,sigma)
    a = zeros(1,100000);
    for i = 1:100000
        u1 = rand;
        u2 = rand;
        r = sqrt(-2 * log(u1));
        v = 2 * pi * u2;
        x = r * cos(v);
        y = r * sin(v);
        a(i) = mu + sigma * x;
    end
end
  
function a = genpoisson(lambda)
    a = zeros(1,10000);
    for i = 1:10000
        u = rand;
        a(i) = -1/lambda * log(1-u);
    end
end

function X = accreja(n)
f = @(x) (100+4/5*x^3-6*x^2)/1000
X = zeros(1,n);
for i=1:n
    accept = false;
    while accept == false
        gx = 0.3 * rand;
        range = 10 * rand;
        if gx <= f(range)
            X(i) = range;
            accept = true;
        end
    end
end
end

function [X,rate] = accrejc(n)
f = @(x) (100+4/5*x^3-6*x^2)/1000;
g = @(x) 0.05 + 0.01*(x-5)^2;
X = zeros(1,n);
rej = 0
for i=1:n
    accept = false;
    while accept == false
        %Hold up, do you have to use Inverse Method?
        x = accrejg();
        gx = g(x) * rand;
        if gx <= f(x)
            X(i) = x;
            accept = true;
        else
            rej = rej + 1;
        end
    end
end
rate = rej/(rej+10000);
end

function X = accrejg()
f = @(x) 0.05 + 0.01*(x-5)^2;
accept = false;
while accept == false
    gx = 0.3 * rand;
    range = 10 * rand;
    if gx <= f(range)
        X = range;
        accept = true;
    end
end
end

function [I,se] = meanvalue(N)
    a = 0;
    b = 1;
    % We sample N points uniformly distributed in [a,b]
    x = rand(1,N) * (b-a) + a;
    % This is the function of which we want to compute the integral
    f = @(x) x.^3 .* sqrt(1-x);
    % And this is the Mean-Value Monte Carlo estimate
    I = (b - a) / N * sum(f(x));
    se = sqrt((b-a)^2/N*var(f(x)));
end

function [I,se] = importance_sampling(N)
    % Sample N points uniformly in [0,1]
    u = rand(N,1);
    % Use the inverse function method to
    % obtain samples from f(x) = 5 * x^4
    x = u .^ (1/5);
    % Compute the Importance Sampling estimate
    I = 1.0/N * sum((x.^3 .* sqrt(1-x)) ./ (5 * x.^4));
    se = sqrt(var((x.^3 .* sqrt(1-x)) ./ (5 * x.^4)) / N)
end

function [I,se] = isf1(N)
    % Sample N points uniformly in [0,1]
    u = rand(N,1);
    % Use the inverse function method to
    % obtain samples from f(x) = 4 * x^3
    x = u .^ (1/4);
    % Compute the Importance Sampling estimate
    I = 1.0/N * sum((x.^3 .* sqrt(1-x)) ./ (4 * x.^3));
    se = sqrt(var((x.^3 .* sqrt(1-x)) ./ (4 * x.^3)) / N)
end

function [I,se] = isftwo(N)
    % Sample N points uniformly in [0,1]
    u = rand(N,1);
    % Use the inverse function method to
    % obtain samples from f(x) = 3 * x^2
    x = u .^ (1/3);
    % Compute the Importance Sampling estimate
    I = 1.0/N * sum((x.^3 .* sqrt(1-x)) ./ (3 * x.^2));
    se = sqrt(var((x.^3 .* sqrt(1-x)) ./ (3 * x.^2)) / N)
end

function [I,se] = isfthree(N)
    % Sample N points uniformly in [0,1]
    u = rand(N,1);
    % Use the inverse function method to
    % obtain samples from f(x) = 2x
    x = u .^ (1/2);
    % Compute the Importance Sampling estimate
    I = 1.0/N * sum((x.^3 .* sqrt(1-x)) ./ (2 * x));
    se = sqrt(var((x.^3 .* sqrt(1-x)) ./ (2 * x)) / N)
end

function x = rw(p, N)
  x = zeros(N, 1);
  x(1) = 0;
  for i = 2:N
    u = rand;
    if u < p
        x(i) = x(i-1) + 1;
    else
        x(i) = x(i-1) - 1;
    end
  end
end

function p = bindist(x,N)
    p = factorial(N)/(factorial((N-x)/2)*factorial((x+N)/2))*0.5^N;
end

function z = pbindist(x,N,p)
    z = factorial(N)/(factorial((N-x)/2)*factorial((x+N)/2))*p^((N+x)/2)*(1-p)^((N-x)/2);
end