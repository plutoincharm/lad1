function y_n_plus1 = fun_explicitEuler(handle_f,y_n,u,dt)

%%% nargin = num of arguments passed to function

    y_n_plus1 = y_n + dt*feval(handle_f,y_n,u);
