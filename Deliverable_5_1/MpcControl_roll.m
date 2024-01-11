classdef MpcControl_roll < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(H/Ts); % Horizon steps
            N = N_segs + 1;      % Last index in 1-based Matlab indexing
            
            [nx, nu] = size(mpc.B);
            
            % Steady-state targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % Horizon and cost matrices
            Q=100*eye(nx);
            R=0.5;
            
            % Slack variables and weight for the quadratic penalty
            epsilon = sdpvar(nx, N);
            S = 100000*eye(nx);

            %Set determination
            x_min = [-Inf; -pi]; 
            x_max = [Inf; pi];
            u_min = [-20]; 
            u_max = [20];
            [~,Qf] = dlqr(mpc.A,mpc.B,Q,R);
            
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
                {X(:,1), x_ref, u_ref}, {U(:,1), X, U});
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state targets
            nx = size(mpc.A, 1);
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            
            % Set determination
            x_min = [-Inf; -pi]; 
            x_max = [Inf; pi];
            u_min = [-20]; 
            u_max = [20];

            % Objective and constraints definition
            obj = 0;
            con = (u_min <= us <= u_max) + ...
                (x_min <= xs <= x_max) + ...
                (xs == mpc.A*xs + mpc.B*us) + ...
                (ref == mpc.C*xs + mpc.D*us);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
