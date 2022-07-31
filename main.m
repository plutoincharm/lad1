clc; clearvars;

global m r L MI g nx ny lam Q R Qf XT xtraj utraj m Nt Denom

L = [1,2];
r = L./2;
m = [1,2]; %[3,2,1]; %[1,2,3]; %
MI = m.*L.^3./12;
g = 9.81;

Q = eye(4);
R = 0.1*eye(2);
Qf = 10*eye(4);
%h = 0.05 ;
Nx = 4;
Nu  = 2;
Tf = 10;
dt = 1;
Nt = round(Tf/dt)+1


th1 = deg2rad(25);
th2 = deg2rad(20);
%%%%%%% var order: (th1,th2,omg1,omg2) rad and rad/s

x0 = [th1;th2;0.2;0.2];
XT=[2*pi;2*pi;0.7;0.7];

xtraj = x0;
utraj = abs(1*randn(Nu,Nt-1));

handle_f = @funDynamic; 
Denomall = [Denom];
%Initial Rollout
for k = 1:(Nt-1)
    %xtraj[:,k+1] .= dynamics_rk4(xtraj[:,k],utraj[k])
    %xtraj(:,k+1)  = funDynamic(xtraj(:,k),utraj(:,k)) 
    xtraj(:,k+1)  = funExplicitEuler(handle_f,xtraj(:,k),utraj(:,k),dt);  % xnd u vector at a k timestamp
    xtrajd(:,k)  = rad2deg(xtraj(:,k+1)); 
       Denomall = [Denomall,Denom];
end
J = funCost(xtraj,utraj,Nt);


%%%%%%%%%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
p = ones(Nx,Nt);
P = zeros(Nx,Nx,Nt);
d = ones(Nu,Nu,Nt-1);
K = zeros(Nu,Nx,Nt-1);
dJ = 1.0;

xn = zeros(Nx,Nt);
un = zeros(Nt-1);

gx = zeros(Nx);
gu = 0.0;
Gxx = zeros(Nx,Nx);
Guu = 0.0;
Gxu = zeros(Nx);
Gux = zeros(Nx);

iter = 0;
while max(abs(d(:))) >  1e-3
    
    iter = iter +  1 ;   

    dJ = backward_pass(p,P,d,K,Qf,XT,xtraj,utraj,Nt,Q,R)          %Backward Pass

    %Forward rollout with line search
    xn(:,1) = xtraj(:,1);
    alpha = 1.0;

    for k = 1:(Nt-1)
        un(k) = utraj(k) - alpha*d(k) - dot(K(:,:,k),xn(:,k)-xtraj(:,k))
        xn(:,k+1) = dynamics_rk4(xn(:,k),un(k))
    end
    Jn = cost(xn,un);
    
    while isnan(Jn) || Jn > (J - 1e-2*alpha*dJ)
        alpha = 0.5*alpha
        for k = 1:(Nt-1)
            un(k) = utraj(k) - alpha*d(k) - dot(K(:,:,k),xn(:,k)-xtraj(:,k))
            xn(:,k+1) = fun_explicitEuler(handle_f,xtraj(:,k),utraj(k),dt);
            %xn(:,k+1) = dynamics_rk4(xn(:,k),un(k)
        end
        Jn = cost(xn,un);
    end
    
    J = Jn;
    xtraj = xn
    utraj = un
   %{ 
    #display(iter)
    #display(α)
    #display(maximum(abs.(d[:])))
   %}
    
end




%%%%%%%%%%%%%%%%   Backward pass      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dJ =  backward_pass(p,P,d,K,Qf,XT,xtraj,utraj,Nt,Q,R)

    dJ = 0.0;
    p(:,Nt) = Qf*(XT- xtraj(:,Nt)); 
    P(:,:,Nt) = Qf;
    
    for k = (Nt-1):-1:1
        
            %Calculate derivatives
            q = Q*(XT-xtraj(:,k));
            r = R*utraj(k);
    
            [A,B] = funJac(xtraj(:,k), utraj(:,k)); % both utraj and xtraj matrix
        
     
            gx = q + A'*p(:,k+1);
            gu = r + B'*p(:,k+1);
    
             %iLQR (Gauss-Newton) version
             Gxx = Q + A'*P(:,:,k+1)*A;
             Guu = R + B'*P(:,:,k+1)*B;
             Gxu = A'*P(:,:,k+1)*B;
             Gux = B'*P(:,:,k+1)*A;     
             
                %{
        beta = 0.1
        while !isposdef(Symmetric([Gxx Gxu; Gux Guu]))
            Gxx += A'*β*I*A
            Guu += B'*β*I*B
            Gxu += A'*β*I*B
            Gux += B'*β*I*A
            β = 2*β
            #display("regularizing G")
            #display(β)
        end
        %}
        
            d(:,:,k) = Guu\gu;
            K(:,:,k) = Guu\Gux;
    
          %  p(:,k) = dropdims(gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(k) - Gxu*d(k), dims=2);
            P(:,:,k) = Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux;
    
            dJ = dJ +  gu'*d(k);

    end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = funCost(xtraj,utraj,Nt)

    J = 0.0;
    for k = 1:(Nt-1)
        J = J + funStagecost(xtraj(:,k),utraj(:,k));
    end
    J = J+ funCost_terminal(xtraj(:,Nt));
  
end








  
  
  