classdef MpcControl_x < MpcControlBase
    
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(H/Ts); % Horizon steps
            N = N_segs + 1;      % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.2)
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
            Q=eye(nx);
            Q(1,1)=8*Q(1,1);
            R=20;
            
            %Calculate the Terminal Set
            sys = LTISystem('A',mpc.A,'B',mpc.B);
            sys.x.min = [-Inf; -0.1745; -Inf; -Inf]; 
            sys.x.max = [Inf; 0.1745; Inf; Inf];
            sys.u.min = [-0.26]; 
            sys.u.max = [0.26];
            sys.x.penalty = QuadFunction(Q); 
            sys.u.penalty = QuadFunction(R);
            Xf = sys.LQRSet;
            Qf = sys.LQRPenalty;
            Ff=Xf.A;
            ff=Xf.b;

            %Display the Terminal Set
            figure;
            subplot(1,3,1)
            grid on;
            Xf.projection(1:2).plot();
            xlabel('angular velocity w_y'); ylabel('angle beta');               
            subplot(1,3,2)
            grid on;
            Xf.projection(2:3).plot();
            xlabel('angle beta'); ylabel('velocity v_x');            
            subplot(1,3,3)
            grid on;
            Xf.projection(3:4).plot();
            xlabel('velocity v_x'); ylabel('position x');
            sgtitle('x Terminal set');

            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            con = (X(:,2) == mpc.A*X(:,1) + mpc.B*U(:,1)) + (sys.u.min <= U(:,1) <= sys.u.max);
            obj = U(:,1)'*R*U(:,1);
            
            for i = 2:N-1
                con = con + (X(:,i+1) == mpc.A*X(:,i) + mpc.B*U(:,i));
                con = con + (sys.x.min <= X(:,i) <= sys.x.max) + (sys.u.min <= U(:,i) <= sys.u.max);
                obj = obj + X(:,i)'*Q*X(:,i) + U(:,i)'*R*U(:,i);
            end

            con = con + (Ff*X(:,N) <= ff);
            obj = obj + X(:,N)'*Qf.weight*X(:,N);

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
            
            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            obj = 0;
            con = [xs == 0, us == 0];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
