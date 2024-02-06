function j = z_axis_mpc(K,dt,p_0,v_0,a_0,pt,vt,at)
 %Implement your code here
    %% Weights for the cost function
    w1 = 5; % Weight for position
    w2 = 1;   % Weight for velocity
    w3 = 1;   % Weight for acceleration
    w4 = 1;   % Weight for jerk
    w5 = 100;
    w6 = 100;
    w7 = 100;
    w8 = 100;
    %% Construct the prediction matrix
    [Tp, Tv, Ta, Bp, Bv, Ba] = getPredictionMatrix(K,dt,p_0,v_0,a_0);
    %% Construct the optimization problem
    H = blkdiag(w4*eye(K) + w1*(Tp'*Tp) + w2*(Tv'*Tv) + w3*(Ta'*Ta),w5*eye(K),w6*eye(K),w7*eye(K),w8*eye(K));
    F = [w1 * (Bp - pt)'*Tp + w2 * (Bv - vt)'*Tv + w3 * (Ba - at)'*Ta zeros(1,K) zeros(1,K) zeros(1,K) zeros(1,K)];

    A = [Tv            -eye(K) zeros(K) zeros(K) zeros(K);%velovity constrains
        -Tv            zeros(K) -eye(K) zeros(K) zeros(K);
         Ta            zeros(K) zeros(K) -eye(K) zeros(K);%acceleration constrain
        -Ta            zeros(K) zeros(K) zeros(K) -eye(K);
        eye(K)         zeros(K) zeros(K) zeros(K) zeros(K);%jerk constrain
       -eye(K)         zeros(K) zeros(K) zeros(K) zeros(K);
       zeros(size(Tv)) -eye(K) zeros(K)  zeros(K) zeros(K);%L1
       zeros(size(Tv)) zeros(K) -eye(K) zeros(K) zeros(K);%L2
       zeros(size(Ta)) zeros(K) zeros(K) -eye(K) zeros(K);%L3
       zeros(size(Ta)) zeros(K) zeros(K) zeros(K) -eye(K);%L4
       ];

    b = [6*ones(K,1)-Bv; %velovity constrains
         1*ones(K,1)+Bv;
         3*ones(K,1)-Ba; %acceleration constrain
         1*ones(K,1)+Ba;
         2*ones(K,1); %jerk constrain
         2*ones(K,1);
         zeros(K,1); %L1
         zeros(K,1); %L2
         zeros(K,1); %L3
         zeros(K,1); %L4
         ];

    %% Solve the optimization problem for z axis
    J = quadprog(H,F,A,b);

    %% Apply the control for z axis (first element of the optimal jerk)
    j = J(1);
end