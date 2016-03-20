function [c1,c2,error]=twiddle(tolerance)
n_params=2;
dparams=[0.1 0.1];
params=[0.5 2];
best_error=Particle_filter2D(100,500,[-0.05;0.001;0.7;-0.055],params(1),params(2));
n=0;
while sum(dparams)>tolerance,
    for i=1:n_params,
        params(i) = params(i) + dparams(i);
        err=Particle_filter2D(100,500,[-0.05;0.001;0.7;-0.055],params(1),params(2));
        if err < best_error,
            best_error=err;
            dparams(i)=dparams(i)*1.1;
        else
            params(i)=params(i)-2.0*dparams(i);
            err=Particle_filter2D(100,500,[-0.05;0.001;0.7;-0.055],params(1),params(2));
            if err<best_error,
                best_error=err;
                dparams(i) = dparams(i) * 1.1;
            else
                params(i)=params(i)+dparams(i);
                dparams(i)=dparams(i)*0.9;
            end
        end
    end
    n=n+1
    params
    dparams
    best_error
end
c1=params(1);
c2=params(2);
error=best_error;