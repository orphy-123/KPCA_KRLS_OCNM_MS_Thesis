
clc
clear all
tic
jpgFiles = dir('C:\MATLAB7\work\*.jpeg');
for k = 1:length(jpgFiles)
  filename = jpgFiles(k).name;
  a= double(imread(filename)); 
   a1=dwt2(a,'haar');
%     a1 = fft2(a);
%   a2=imresize(abs(a1),0.1,'bil');
 a2=imresize(a1,0.1,'bil');

[m n]=size(a2);
a2_expand=[];
for r=1:m
a2_expand=[a2_expand a2(r,:)];
end
X(k,:)=[a2_expand];
end

[T f] = size(X);


%KRLS online anomaly detection algorithm.
%Requires the following inputs:
% Mandatory: data matrix (X), lower threshold (nu1), upper threshold (nu2).
% Parameters for dropping obsolete elements; default is d = 0.9, L = 100.
% Parameters for resolving orange alarm; default is epsilon = 0.2, el = 20.
% Parameters for resetting P; default is R = 10000, r = 1.
% Forgetting factor; default is gamma = 1;
%Yields the following outputs:
% Always: Red1 and Red2 alarm positions.
% If desired: records of delta and prediction error.

%function [Red1_out, Red2_out deltaStore_out Error_out] = KRLS(X, nu1, nu2, d, L, epsilon, el, R, r, gamma)

%if nargin < 11 = 1; end %Needed only if using Gaussian kernel
gamma = 1; %Forgetting factor
r = 1;  %Parameters for resetting P
R = 10000; 
el = 20;  %Parameters for resolving orange alarm
epsilon = 0.2; 
L = 100;  %Parameters for dropping obsolete elements
d = 0.9; 

%X = P; %May be X=B or X=F to refer to Bytes or Flows, instead of Packets
X = X./repmat(sqrt(sum(X.*X,2)+eps),1,size(X,2)); %normalize to unit circle (i.e. divide by norm)
Y = sum(X,2); %Add after normalizing
[T f] = size(X);

nu1 = 0.001; nu2 = 0.0012; %Threshholds
%d = 0.9; L = 100; %Parameters for dropping obsolete elements
%epsilon = 0.20; el = 20; %Parameters for resolving orange alarm
%R = 10000; r = 1; %Parameters for resetting P
%gamma = 1; %Forgetting factor
%sigma = 1; %Needed only if using Gaussian kernel

Red1 = []; Red2 = [];%Clear alarms
Orange = []; x_Orange = []; %Store x in timesteps when Orange alarm is raised

% Initialize %
t = 1;
x = X(t,:)';
y = Y(t);
k11 = kernel(x, x);
K_tilde = [k11];
K_tilde_inv = [1/k11];
Dictionary = [x];
index_m = [t]; %Keeps track of timesteps when elements are added (+2), deleted (-1) or no change to D (0); for debugging only
Orange = [Orange t];
x_Orange = [x_Orange x];
drop_index = [0];
P=[1]; %P=inv(A'A)
m=1;
m_t(t) = m; %Keep track of m, for debugging only
index_m(t) = 2; %index_m(t)=2 implies x(t) is being added to Dictionary
alpha = y(t)/k11;
deltaStore(1) = nu1+eps;  %For debugging

% % Evaluate y_hat %
y_hat = zeros(T,1);
Error = zeros(1,T);
for j=1:m
    y_hat(t) = y_hat(t) + alpha(j)*kernel(Dictionary(:,j),x);
end %for j=1:m
Error(t) = (Y(t)-y_hat(t))/Y(t)*100;

clear Lambda lambda dotProd;
Lambda = kernel(Dictionary(:,j),x);

%Keep track of all dot product (kernel) values; for debugging only
for j=1:m
    dotProd(t,j) = kernel(Dictionary(:,j),x);
end %for j=1:m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:T
    %t=t+1;
    x = X(t,:)';
    y = Y(t);

    % Evaluate current measurement %
    k_tilde = zeros(m,1);
    for j=1:m
        k_tilde(j) = kernel(Dictionary(:,j),x); %Computing k_tilde{t-1}
    end %for j=1:m
    a = K_tilde_inv*k_tilde;
    delta = kernel(x,x) - k_tilde'*a;
%    deltaCheck = a'*K_tilde*a - 2*a'*k_tilde + kernel(x,x); %Verify delta; not part of algorithm
    deltaStore(t) = delta;  %Keep track of delta, for debugging only
    if t>L
        Lambda = [Lambda(2:end,1:end) ; ceil(k_tilde'-repmat(d,1,m))]; %Append with 1 or 0
    else
        Lambda = [Lambda ; ceil(k_tilde'-repmat(d,1,m))]; %Append with 1 or 0
    end %if t>L

    if (delta>=nu1 & delta<nu2) %Orange alarm, add to Dictionary
        x_Orange = [x_Orange x];
        Orange = [Orange t];
        Dictionary = [Dictionary x];
        drop_index = [drop_index 0];
        a_tilde = a;
        K_tilde_inv = [ (delta*K_tilde_inv+a_tilde*a_tilde') (-1*a_tilde) ; (-1*a_tilde') (1) ] / delta;
        K_tilde = [ K_tilde k_tilde ; k_tilde' kernel(x,x) ];
        if t>L
            lambda = [zeros(L-1,1) ; 1];
        else
            lambda = [zeros(t-1,1) ; 1];
        end %if t>L
        Lambda = [Lambda lambda];

        a = [zeros(m-1,1) ; 1];
        P = [ P zeros(m,1) ; zeros(m,1)' gamma ]/gamma;
        alpha = [ (gamma^(-0.5)*alpha - a_tilde*(y-gamma^(-0.5)*k_tilde'*alpha)/delta)  ;  ((y-gamma^(-0.5)*k_tilde'*alpha)/delta) ];
        m=m+1;
        m_t(t) = m;
        index_m(t) = 2; %Element added to D in this timestep
    else %delta<nu1 or delta>=nu2, Dictionary unchanged
        if delta>nu2 %Red1 alarm
            Red1 = [Red1 t];
        end %if delta>nu2
        %K_tilde = K_tilde;
        %K_tilde_inv = K_tilde_inv;
        q = (P*a) / (gamma+a'*P*a);
        P = (1/gamma)*[P - q*a'*P];
        alpha = alpha + K_tilde_inv*q*(y-k_tilde'*alpha);
        %m = m;
        m_t(t) = m;
        index_m(t) = 0; %No change to D in this timestep
    end %delta > nu1 & delta < nu2 


    %Keep track of all dot product (kernel) values; for debugging only
    for j=1:m
        dotProd(t,j) = kernel(Dictionary(:,j),x);
    end %for j=1:m
    
    % Process previous orange alarm %
    if t>el & sum(Orange==t-el)==1 %means orange alarm at timestep t-el
        %Identify Dictionary element j corr. to the orange alarm at timestep t-el
        for j=1:m
            if x_Orange(:,Orange==t-el)==Dictionary(:,j)
                break;
            end %if x_Orange(:,Orange==t-el)==Dictionary(:,j)
        end %for j=1:m

        if sum(Lambda(end-el+1:end,j)) <= epsilon*el
            %Orange turns Red
            Red2 = [Red2  Orange(Orange==t-el)]; %Red2 alarm
            x_Orange(:,Orange==t-el) = [];
            Orange(Orange==t-el) = [];
            drop_index = [zeros(1,j-1) 1 zeros(1,m-j) ];
        else
            %Orange turns green
            x_Orange(:,Orange==t-el) = [];
            Orange(Orange==t-el) = [];
        end %if size(find(Lambda(end-el+1:end,j)<d),1) >= 0.80*25
    end %if t>el & sum(Orange==t-el)==1


    % Remove obsolete elements %
    for j=1:m
        %Dropping condition: kernel exists for past L timesteps, and is always < d
        if  ( t>L & sum(Lambda(1:end,j))==0 )
            drop_index(j) = 1;
        end %if  ( t>L & gt(Lambda(:,j),0) & lt(Lambda(:,j),d) )
    end %for j=1:m

    % DropElement(p) %
    if ( find(drop_index==1) & m>1 & t>r )
        t;
        p = min(find(drop_index==1)); %Drop Dictionary element # p
        %Reorganize K_tilde_p and K_tilde_inv_p, with p'th row/col moved to the end
        K_tilde = [ K_tilde(1:p-1,1:p-1) K_tilde(1:p-1,p+1:m) K_tilde(1:p-1,p)  ;  K_tilde(p+1:m,1:p-1) K_tilde(p+1:m,p+1:m) K_tilde(p+1:m,p)  ;  K_tilde(p,1:p-1) K_tilde(p,p+1:m) K_tilde(p,p) ];
        K_tilde_inv = [ K_tilde_inv(1:p-1,1:p-1) K_tilde_inv(1:p-1,p+1:m) K_tilde_inv(1:p-1,p)  ;  K_tilde_inv(p+1:m,1:p-1) K_tilde_inv(p+1:m,p+1:m) K_tilde_inv(p+1:m,p)  ;  K_tilde_inv(p,1:p-1) K_tilde_inv(p,p+1:m) K_tilde_inv(p,p) ];
        delta_p = 1/(K_tilde_inv(m,m));
        a_tilde_p = -delta_p*[K_tilde_inv(1:m-1,m)];
        K_tilde_inv = K_tilde_inv(1:m-1,1:m-1)-a_tilde_p*a_tilde_p'/delta_p;
        alpha = alpha - (1/delta_p)*[a_tilde_p*a_tilde_p'  -a_tilde_p  ;  -a_tilde_p'  1] *K_tilde*alpha;
        alpha = alpha(1:m-1);
        K_tilde = K_tilde(1:m-1,1:m-1);
        Dictionary(:,p) = [];
        drop_index(p) = [];
        Lambda(:,p) = [];
        dotProd(:,p) = []; %%%
        m=m-1;
        m_t(t) = m;
        index_m(t) = -1; %Element deleted from D in this timestep

        % Reset P %
        P = R*eye(m);
        for i_r=1:r
            k_tilde = zeros(m,1);
            for j=1:m
                k_tilde(j) = kernel(Dictionary(:,j),X(t-i_r,:)'); %Computing k_tilde{t-1}
            end %for j=1:m
            a = K_tilde_inv*k_tilde;
            q = (P*a) / (gamma+a'*P*a);
            P = (1/gamma)*[P - q*a'*P];
            alpha = alpha + K_tilde_inv*q*(Y(t-i_r)-k_tilde'*alpha);
        end %for i_r=1:r-1
     end %if ( find(drop_index==1) & m>1 & t>r )

%     % Evaluate y_hat %
%     for j=1:m
%         y_hat(t) = y_hat(t) + alpha(j)*kernel(Dictionary(:,j),x);
%     end %for i=1:m
%     Error(t) = (Y(t)-y_hat(t))/Y(t)*100;

end %for t=2:T


    Red1_out = Red1; Red2_out = Red2; 

    Red1_out = Red1; Red2_out = Red2; deltaStore_out = deltaStore;

    Red1_out = Red1; Red2_out = Red2; deltaStore_out = deltaStore; Error_out = Error;

figure(1)
stem(deltaStore_out)

for b=2:length(jpgFiles)
    if deltaStore_out(b)>1*10^-3
        figure(b)
        imshow(jpgFiles(b).name)
    end
end

toc