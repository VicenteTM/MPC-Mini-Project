    classdef MpcControl_zwe < MpcControlBase
    properties
        A_bar, B_bar, C_bar % Augmented system for disturbance rejection
        L                   % Estimator gain for disturbance rejection
    end
    
    methods
        function mpc = MpcControl_zwe(sys, Ts, H)
            mpc = mpc@MpcControlBase(sys, Ts, H);
            
            [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
        end
        
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   d_est        - disturbance estimate
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(H/Ts); % Horizon steps
            N = N_segs + 1;      % Last index in 1-based Matlab indexing
            
            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.3)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar(1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % Horizon and cost matrices
            Q = 100*eye(nx);
            Q(2,2)=2*Q(2,2);
            R = 3;
            [~,Qf] = dlqr(mpc.A,mpc.B,Q,R);
            
            % Slack variables and weight for the quadratic penalty
            epsilon = sdpvar(nx, N);
            S = 100000*eye(nx);

            %Set definition
            x_min = [-Inf; -Inf]; 
            x_max = [Inf; Inf];
            u_min = [50-56.6667]; 
            u_max = [80-56.6667];

            % Objective and constraints definition
            con = ((X(:,2) - x_ref) == mpc.A*(X(:,1) - x_ref) + mpc.B*(U(:,1) - u_ref)) + ...
                ((x_min - epsilon(:,1)) <= (X(:,1) - x_ref) <= (x_max + epsilon(:,1))) + ...
                (u_min <= (U(:,1) - u_ref) <= u_max) +  (epsilon(:,1) >= zeros(nx,1));
            obj = (U(:,1) - u_ref)'*R*(U(:,1) - u_ref) + (epsilon(:,1)'*S*epsilon(:,1));
            
            for i = 2:N-1
                con = con + ((X(:,i+1) - x_ref) == mpc.A*(X(:,i) - x_ref) + mpc.B*(U(:,i) - u_ref));
                con = con + ((x_min - epsilon(:,i)) <= (X(:,i) - x_ref) <= (x_max + epsilon(:,i))) + ...
                    (u_min <= (U(:,i) - u_ref) <= u_max) + (epsilon(:,i) >= zeros(nx,1));
                obj = obj + (X(:,i) - x_ref)'*Q*(X(:,i) - x_ref) + ...
                    (U(:,i) - u_ref)'*R*(U(:,i) - u_ref) + (epsilon(:,i)'*S*epsilon(:,i));
            end

            obj = obj + (X(:,N)-x_ref)'*Qf*(X(:,N)-x_ref) + (epsilon(:,N)'*S*epsilon(:,N));
            con = con + ((x_min - epsilon(:,N)) <= (X(:,N)-x_ref) <= (x_max + epsilon(:,N))) + ...
                (epsilon(:,N) >= zeros(nx,1));

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref, d_est}, {U(:,1), X, U});
        end
        
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);
            
            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.3)
            ref = sdpvar;
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            
            % Defining the sets
            u_min = 50-56.6667;
            u_max = 80-56.6667;
            
            Bd = [1;
                  0];
              
            Cd = 1;
            
            % Defining the objective and the constraints
            obj = 0;
            con = (u_min <= us <= u_max) + ...
                (xs == mpc.A*xs + mpc.B*us + Bd* d_est) + ...
                (ref == mpc.C*xs + Cd*d_est);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
        end
        
        
        % Compute augmented system and estimator gain for input disturbance rejection
        function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
            
            %%% Design the matrices A_bar, B_bar, L, and C_bar
            %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
            %%% converges to the correct state and constant input disturbance
            %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            nx   = size(mpc.A,1);
            nu   = size(mpc.B,2);
            ny   = size(mpc.C,1);
            
            Bd = [1;
                  0];
              
            Cd = 1;
            
            A_bar = [mpc.A, Bd;       
                     zeros(1,nx),1];
            B_bar = [mpc.B;zeros(1,nu)];
            C_bar = [mpc.C,Cd];        
    
            L = -place(A_bar',C_bar',[0.3, 0.4, 0.5])';   %Maybe would need to change poles
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
    end
end
