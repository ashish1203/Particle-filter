function[g,b] = re(xp,pdf) % resampling
    c(1)=0;
    for n=2:500
        c(n)=c(n-1)+pdf(n);
    end 
    k=1;
    u=zeros(500,1);
    u=1/500*rand(500,1);
    
    for j=1:500
        
        u(j)=u(1)+(j-1)/500;
        while u(j)>c(k)
            k=k+1;
            if k>500
                k=k-1;
                break;
            end 
        end
    
        g(:,j)=xp(:,k);
        b(j)=1/500;
     %   p(:,496:500)=xp(:,500);
      %  b(496:500)=1/500;
    end 
       

    
	