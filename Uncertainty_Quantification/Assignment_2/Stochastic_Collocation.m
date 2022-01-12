% this script integrates the SEIR model M times
% with tau distributed according to Clenshaw-Curtis nodes.

figure(1); clf;

% ls will be the list of Q-values.
ls = [];

M = 14;
% Curtis-clenshaw nodes (cos(i*pi/M)):
tt = linspace(0,M,M);
% Rescale to [1,14]
x = rescale(cos(tt*pi/M),1,14);

for n=1:M

    R0=2.2;
    T=9; 
    %tau=1.0+13*betarnd(2,2,1,1);
    tau = x(n);
    %SEIRmodel is the given file (unchanged)
    Xoutput = SEIRmodel(R0,T,tau);
    
    Q = max(Xoutput(:,3)); % maximum of timeseries of I
    ls = [ls, Q];

% Visualization of the nodes according to maxima of I.
% Also uncomment the top line in order to see figure 1.
    figure(1);% subplot(2,1,1);
%     plot(Xoutput(:,3),'-r'); hold on; % plot timeseries of I
%     xlabel('timesteps'); ylabel('I')
%     subplot(2,1,2);
    plot(tau,Q,'or');
    plot(x,-0.5*ones(size(x)),'xk');
    hold on;
    xlabel('\tau'); ylabel('Q')
    
end
%%
% Create grid from 1 to 14, with steps of 0.1
% for interpolation. Can increase to smoothen curve.
dz=0.1;
z=(1:dz:14);
J=length(z)-1;

% fa will be approximated function
fa=zeros(size(z));

for j=1:length(z)
for i=1:M
    qi=lagrangebasisfunc(z(j),x,i);
    fa(j)=fa(j)+ls(i)*qi;
end
end

fa_mean = mean(fa);
fa_std = std(fa);

% create figure for approximated function, with interpolated nodes
% and their x-axis value visible.
figure(2); clf;
plot(z,fa,'-r',x,ls,'or',x,-0.5*ones(size(x)),'xk','LineWidth',2)
legend('approximation','interpolation points','Location','NW')
xlabel('\tau'); ylabel('Q')
title(['Interpolation, M=',int2str(M)])
hold on

%% 
% Monte Carlo Sampling of the SEIR model.

n = 1000;
Q_seir = zeros(1, n);
tau_c = zeros(1, n);
R0 = 2.2;
T = 9;

for i=1:n

    tau = 1+rand*13;

    tau_c(i) = tau;
    
    xoutput = SEIRmodel(R0,T,tau);
    
    Q_seir(i) = max(xoutput(:,3));
    
end

figure
plot(tau_c, Q_seir, '.')

Q_mean = mean(Q_seir);
Q_std = std(Q_seir);

histogram(fa, 20);
xlabel('Q'); ylabel('Numer of nodes')

%%
function y=lagrangebasisfunc(x,xi,i)

% evaluate i-th Lagrange basis function at x, given nodes xi
% y = l_i(x) = \product(j \neq i) (x-x_j)/(x_i-x_j)
% according to the definition in Xiu

y=1;
 
for j=1:length(xi)
    if i~=j
    y = y*(x-xi(j))/(xi(i)-xi(j));
    end
end

end



