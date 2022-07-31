function f = funDynamic(x,u)  % x and u are vector at that time stamp

[A,B] = funJac(x,u);
f  = A*x + B*u;

