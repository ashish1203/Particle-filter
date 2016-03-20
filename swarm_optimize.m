function [xpp,weights]=swarm_optimize(p,q,y,c1,c2)
N=length(q);
bearing_cov=0.01;
V=zeros(4,N);
yp=zeros(1,N);
pbest_weight=q;
pbest_particles=p;
[~,idx]=max(pbest_weight);
gbest=pbest_particles(:,idx);
%repel=zeros(size(p));
% pause;
%  plot(p(1,:),p(3,:),'m*');
%  hold on;
%  pause;
for k=1:5
%     for b=1:N,
%     particle=p(:,b);
%     particle=particle*ones(1,size(p,2));
%     temp=p-particle;
%     temp(2,:)=temp(2,:).*0;
%     temp(4,:)=temp(4,:).*0;
%     temp=temp.^2;
%     dist=sum(temp);
%     variable=dist<0.01;
%     variable(b)=0;
%     temp=(variable'*ones(size(variable,1),size(p,1)))';
%     variable=p.*temp;
%     repel=repel+(variable-particle).*temp;
%     end
    personal= (c1*rand(size(p)).*(pbest_particles-p));
    globals = (c2*rand(size(p)).*((gbest*ones(1,N))-p));
% repel = (c3*rand(size(p))).*repel;
    V = V + personal +  globals ;%+ repel;
    %repel=zeros(size(p));
    p = p + V;
%     plot(p(1,:),p(3,:),'ro');
%     hold on;
%     plot(gbest(1),gbest(3),'gs');
%     pause;
%     hold off;
    for i=1:N
        yp(i)=h(p(:,i),v);
    end
    d=-(yp-y);
    q = 1/sqrt(2*pi*bearing_cov)*exp(-(d.^2)/(2*bearing_cov));
    q=q./sum(q);
    grt_logic=q>pbest_weight;
    pbest_weight=(pbest_weight-pbest_weight.*(grt_logic))+(q.*(grt_logic));
    grt_logic=((grt_logic')*ones(size(grt_logic,1),size(p,1)))';
    pbest_particles=pbest_particles.*(~(grt_logic))+p.*(grt_logic);
    pbest_weight=pbest_weight./sum(pbest_weight);
    [~,idx]=max(pbest_weight);
    gbest=pbest_particles(:,idx);
end

xpp=pbest_particles;
weights=pbest_weight;
% xp(:,1:50)=gbest*ones(1,50);
% temp2=sort(pbest_weight,'descend');
% weights(1:50)=pbest_weight(1);
% p3=51;
% 
% for p1=2:451
%     for p2=1:N
%         if p3<=N
%         if temp2(1,p1)==pbest_weight(1,p2)
%             xp(:,p3)=pbest_particles(1,p2);
%             weights(p3)=pbest_weight(p2);
%             p3=p3+1;
%         end
%         end
%     end
% end


% for p=1:5
%     for r=1:500
%         V(1,r)=V(1,r)+c1*rand*(pbest_particles(1,r)-X_particles(j+1,r))+c2*rand*(gbest-X_particles(j+1,r));
%         X_particles(j+1,r)=X_particles(j+1,r)+V(1,r);
%         y=(X_particles(j+1,r)^2)/20;
%         weight(1,r)=exp(-0.5*(M_N^(-1))*((Z(j+1,1)-y)^2));
%         if weight(1,r)>pbest_weight(1,r)
%             pbest_weight(1,r)=weight(1,r);
%             pbest_particles(1,r)=X_particles(1,r);
%         end
%     end
%     temp=max(weight);
%       j1=1;
%       while weight(1,j1)~=temp
%           j1=j1+1;
%       end
%       gbest=pbest_particles(1,j1);
% end
