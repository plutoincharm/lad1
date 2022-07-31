function dJ = funBWDpass(p,P,d,K)
    
global m r L MI g nx ny lam Q R Qf XT x xtraj utraj

    dJ = 0.0;
    p(:,Nt) = Qf*(XT- xtraj(:,Nt)); 
    P(:,:,Nt) = Qf;
    
    for k = (Nt-1):-1:1
        
        %Calculate derivatives
        q = Q*(xtraj[:,k]-xgoal)
        r = R*utraj[k]
    
        [A,B] = funJac(xtraj(:,k), utraj(:,k)); % both utraj and xtraj matrix
        
     
        gx = q + A'*p(:,k+1);
        gu = r + B'*p(:,k+1);
    
        %iLQR (Gauss-Newton) version
        Gxx = Q + A'*P(:,:,k+1)*A
        Guu = R + B'*P(:,:,k+1)*B
        Gxu = A'*P(:,:,k+1)*B
        Gux = B'*P(:,:,k+1)*A
        
     
        
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
        
        d(k) = Guu\gu;
        K(:,:,k) = Guu\Gux;
    
        p(:,k) = dropdims(gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(k) - Gxu*d(k), dims=2)
        P(:,:,k) = Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux
    
        dJ = dJ +  gu'*d(k)
    end
    
    % return dJ
