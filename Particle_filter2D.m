function rmse=Particle_filter2D(T,x,N)
clc;
clear all;
init=[-0.05;0.001;0.7;-0.055];
x=init;
T=30;
N=500;
init=x;
value=zeros(100,T+1);
for loop=1:1
x=init;
X=zeros(size(x,1),T+1);% particles with states
xp=zeros(size(x,1),N);% values for state for each particle
yp=zeros(1,N);
q=ones(1,N)./N;  % Initial weights assigned ...
bearing_cov=0.01;
sigma=[0.5;0.005;0.3;0.01];
p= x*ones(1,N) + ((sigma.^2)*ones(1,N)).*randn(4,N); 
% plot(x(1),x(3),'o');
% hold on;
% plot(p(1,:),p(3,:),'m*');
% axis([-0.25,0.3,-.75,.8]);
% pause(0.2);
X_hat=zeros(size(x,1),T+1);
X_hat(:,1)=(q*p')';  % weighted mean of state....
for t=1:10
    X(:,t)=x;
    x=f(x,[0.0;0.0],0);
    y=h(x,0);
  
    for i=1:N 
        xp(:,i)=f(p(:,i),w,1);
        yp(i)=h(xp(:,i),v);
    end
    p=xp;
    X_hat(:,t+1)=(q*p')';  % weighted particle for predicted.....
%     plot(X(1,1:t),X(3,1:t),'o');
%     hold on;
%     plot(X_hat(1,1:t),X_hat(3,1:t),'gs');
%     plot(p(1,:),p(3,:),'m*');
%     hold off;
%     pause(0.5);
    d=-(yp-y);%
    q = 1/sqrt(2*pi*bearing_cov)*exp(-d.^2/(2*bearing_cov));
    q=q./sum(q);
   % p=re_sample2d(xp,q); % resampling ....
   [g,b]=re(xp,q);
%    xpp=g;
pmodify=1;
pmutate=0.5;
Keep = 4; % elitism parameter: how many of the best habitats to keep from one generation to the next
lambdaLower = 0.0; % lower bound for immigration probabilty per gene
lambdaUpper = 1; % upper bound for immigration probabilty per gene
dt = 1; % step size used for numerical integration of probabilities
I = 1; % max immigration rate for each island
E = 1; % max emigration rate, for each island
P=500;
MinParValue=-21.51;
MaxParValue=24.57;
temp3=zeros(4,500);
particle_keep=zeros(4,2);
weight_Keep=zeros(1,2);
for n = 1 : 500
    Prob(n) = 1 / 500; 
end

                                       %% BBO %%  
   % Initialize the species count probability of each habitat
% Later we might want to initialize probabilities based on cost
% Begin the optimization loop
weight=b;
X_particles=g;
[X_particles,weight] = popsort1(X_particles, weight);
    % Save the best habitats in a temporary array.
    for k = 1 : Keep
        particle_keep(:,k)= X_particles(:,k);
        weight_Keep(1,k) = weight(1,k);
    end
[SpeciesCount]=GetSpeciesCounts(weight, inf);
  % Compute immigration rate and emigration rate for each species count.
    % lambda(i) is the immigration rate for habitat i.
    % mu(i) is the emigration rate for habitat i.
[lambda, mu] = GetLambdaMu(SpeciesCount, I, E);
 
        % Compute the time derivative of Prob(i) for each habitat i.
        for n = 1 : 500
            % Compute lambda for one less than the species count of habitat i.
            lambdaMinus = I * (1 - (SpeciesCount(1,n) - 1) / P);
            % Compute mu for one more than the species count of habitat i.
            muPlus = E * (SpeciesCount(1,n) + 1) / P;
            % Compute Prob for one less than and one more than the species count of habitat i.
            % Note that species counts are arranged in an order opposite to that presented in
            % MacArthur and Wilson's book - that is, the most fit
            % habitat has index 1, which has the highest species count.
            if n < P
                ProbMinus = Prob(n+1);
            else
                ProbMinus = 0;
            end
            if n > 1
                ProbPlus = Prob(n-1);
            else
                ProbPlus = 0;
            end
            ProbDot(n) = -(lambda(n) + mu(n)) * Prob(n) + lambdaMinus * ProbMinus + muPlus * ProbPlus;
        end
        % Compute the new probabilities for each species count.
        Prob = Prob + ProbDot * dt;
        Prob = max(Prob, 0);
        Prob = Prob / sum(Prob); 
    
    % Now use lambda and mu to decide how much information to share between habitats.
    lambdaMin = min(lambda);
    lambdaMax = max(lambda);
    for k = 1 : P
        if rand > pmodify
            continue;
        end
        % Normalize the immigration rate.
        lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda(1,k) - lambdaMin) / (lambdaMax - lambdaMin);
        % Probabilistically input new information into habitat i
% for j = 1 : OPTIONS.numVar
            if rand < lambdaScale
                % Pick a habitat from which to obtain a feature
                RandomNum = rand * sum(mu);
                Select = mu(1);
                SelectIndex = 1;
                while (RandomNum > Select) && (SelectIndex < P)
                    SelectIndex = SelectIndex + 1;
                    Select = Select + mu(SelectIndex);
                end
                temp3(:,k) = X_particles(:,SelectIndex);
            else
                temp3(:,k) = X_particles(:,k);
            end
        
    
    
        % Mutation
        Pmax = max(Prob);
        MutationRate = pmutate * (1 - Prob / Pmax);
        % Mutate only the worst half of the solutions
               %acc. to weight
        for k = round(P/2) : P
            
                if MutationRate(k) > rand
                    temp3(1,k) = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand);
                    temp3(3,k) = floor(MinParValue + (MaxParValue - MinParValue + 1) * rand);
                end
            
        end
    end
    % Replace the habitats with their new versions.
    for k = 1 : P
        X_particles(:,k) = temp3(:,k);
        %X_particles(3,k) = temp3(3,k);
    end

[X_particles,weight] = popsort1(X_particles, weight);
p=X_particles;
q=weight; 
%      figure(1);
%      %imtool(xp(1,:),xp(3,:),'gs');
%      plot(xp(1,:),xp(3,:),'gs');
%      hold off;
%      pause(1);
%      figure(2);
%      plot(X_particles(1,:),X_particles(3,:),'m*');
%      hold off;
%      pause(1);

end
% X_output(j,1)=X_particles(j,1);
%X_output(j,1)=mean(X_particles(j,1:50));
temp=((X(1,:)-X_hat(1,:)).^2 + (X(3,:)-X_hat(3,:)).^2 );
value(loop,:)=temp;
end
temp=sqrt(sum(value)./1);
rmse=(sum(temp)/length(temp));
hold on;
plot(X(1,:),X(3,:),'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
hold on;
