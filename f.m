function x=f(x,w,k)
phi=[1 1 0 0;
     0 1 0 0;
     0 0 1 1;
     0 0 0 1];
tau=[0.5 0;
     1 0;
     0 0.5;
     0 1];
 acov=0.0000005;
 amean=0.0005 + k*sqrt(acov)*randn;
x=phi*x+tau*w+[0;amean;0;0.0];