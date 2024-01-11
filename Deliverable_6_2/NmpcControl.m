classdef NmpcControl < handle
    
    properties
        solver
        nx, nu, N
        nlp_x0
        nlp_lbx, nlp_ubx
        nlp_lbg, nlp_ubg
        nlp_p
        
        T_opt
        sol
        idx
        
        % Delay compensation
        rocket
        expected_delay
        mem_u
        
        % Warmstart
        nlp_lam_x0
        nlp_lam_g0
    end
    
    methods
        function obj = NmpcControl(rocket, tf, expected_delay)
            
            if nargin < 3, expected_delay = 0; end
           
            import casadi.*
            
            N_segs = ceil(tf/rocket.Ts); % MPC horizon
            nx = 12; % Number of states
            nu = 4;  % Number of inputs
            
            % Decision variables (symbolic)
            N = N_segs + 1; % Index of last point
            X_sym = SX.sym('X_sym', nx, N); % state trajectory
            U_sym = SX.sym('U_sym', nu, N-1); % control trajectory)
            
            % Parameters (symbolic)
            x0_sym  = SX.sym('x0_sym', nx, 1);  % initial state
            ref_sym = SX.sym('ref_sym', 4, 1);  % target position
            
            % Default state and input constraints
            ubx = 1e20*ones(nx, 1);
            lbx = -1e20*ones(nx, 1);
            ubu = 1e20*ones(nu, 1);
            lbu = -1e20*ones(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            cost = 0;
            eq_constr = [ ; ];
            ineq_constr = [ ; ];


            % For box constraints on state and input, overwrite entries of
            % lbx, ubx, lbu, ubu defined above

            % |beta| < 75 deg
            lbx(5) = deg2rad(-75);
            ubx(5) = deg2rad(75);

            % |delta1| < 15 deg
            lbu(1) = deg2rad(-15);
            ubu(1) = deg2rad(15);

            % |delta2| < 15 deg
            lbu(2) = deg2rad(-15);
            ubu(2) = deg2rad(15);

            % Pavg in [50,80]
            lbu(3) = 50;
            ubu(3) = 80;

            % |Pdiff| < 20
            lbu(4) = -20;
            ubu(4) = 20;

            % Define the cost function
            Q = diag([200, 200, 100,...
                      1000, 1000, 100, ...
                      100, 100, 100, ...
                      150, 150, 200]); % State cost
            R = diag([0.01, 0.01, 0.01, 0.5]); % Input cost
            
  
            ref = [0, 0, 0, 0, 0, ref_sym(4), 0, 0, 0, ref_sym(1:3)']';
            
            eq_constr = [eq_constr; X_sym(:,1) - x0_sym];

            cost_f = @(x,u) (x-ref)'*Q*(x-ref) + u'*R*u;
            rk4 = @(x,k1,k2,k3,k4) x + rocket.Ts/6*(k1+2*k2+2*k3+k4);
            
            for i = 1:N-1
                
                cost = cost + cost_f(X_sym(:,i), U_sym(:,i));
            
                k1 = rocket.f(X_sym(:,i), U_sym(:,i));
                k2 = rocket.f(X_sym(:,i)+rocket.Ts/2*k1, U_sym(:,i));
                k3 = rocket.f(X_sym(:,i)+rocket.Ts/2*k2, U_sym(:,i));
                k4 = rocket.f(X_sym(:,i)+rocket.Ts*k3,   U_sym(:,i));


                X_p1 = rk4(X_sym(:,i), k1, k2, k3, k4);
                eq_constr = [eq_constr; X_sym(:,i+1) - X_p1];
            end

            [xs, us] = rocket.trim();
            sys = rocket.linearize(xs, us);
            sys = c2d(sys,rocket.Ts);
            [~,Qf] = dlqr(sys.A,sys.B,Q,R); 
            
            cost = cost + (X_sym(:,N)-ref)'*Qf*(X_sym(:,N)-ref); 
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ---- Assemble NLP ------
            nlp_x = [X_sym(:); U_sym(:)];
            nlp_p = [x0_sym; ref_sym];
            nlp_f = cost;
            nlp_g = [eq_constr; ineq_constr];
            
            nlp = struct('x', nlp_x, 'p', nlp_p, 'f', nlp_f, 'g', nlp_g);
            
            % ---- Setup solver ------
            opts = struct('ipopt', struct('print_level', 0), 'print_time', false);
            obj.solver = nlpsol('solver', 'ipopt', nlp, opts);
            
            % ---- Assemble NLP bounds ----
            obj.nlp_x0  = zeros(size(nlp_x));
            
            obj.nlp_ubx = [repmat(ubx, N, 1); repmat(ubu, (N-1), 1)];
            obj.nlp_lbx = [repmat(lbx, N, 1); repmat(lbu, (N-1), 1)];
            
            obj.nlp_ubg = [zeros(size(eq_constr)); zeros(size(ineq_constr))];
            obj.nlp_lbg = [zeros(size(eq_constr)); -inf(size(ineq_constr))];
            
            obj.nlp_p = [zeros(size(x0_sym)); zeros(size(ref_sym))];
            
            obj.nlp_lam_x0 = [];
            obj.nlp_lam_g0 = [];
            
            obj.nx = nx;
            obj.nu = nu;
            obj.N = N;
            obj.T_opt = linspace(0, N * rocket.Ts, N);
            
            obj.idx.X = [1, obj.N * obj.nx];
            obj.idx.U = obj.idx.X(2) + [1, (obj.N-1) * obj.nu];
            obj.idx.u0 = obj.idx.U(1) + [0, obj.nu-1];
            
            % Members for delay compensation
            obj.rocket = rocket;
            obj.expected_delay = expected_delay;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            u_init = [0, 0, 50, 0]';
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.mem_u = repmat(u_init, 1, expected_delay);
        end
        
        function [u, T_opt, X_opt, U_opt] = get_u(obj, x0, ref)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % Delay compensation: Predict x0 delay timesteps later.
            % Simulate x_ for 'delay' timesteps
            x_ = x0;
            rk4 = @(x,k1,k2,k3,k4) x + obj.rocket.Ts/6*(k1+2*k2+2*k3+k4);

            for d = 1:obj.expected_delay
                u_d = obj.mem_u(:,d);

                k1 = obj.rocket.f(x_, u_d);
                k2 = obj.rocket.f(x_+obj.rocket.Ts/2*k1, u_d);
                k3 = obj.rocket.f(x_+obj.rocket.Ts/2*k2, u_d);
                k4 = obj.rocket.f(x_+obj.rocket.Ts*k3,   u_d);
                x_ = rk4(x_, k1, k2, k3, k4);
            end

       
            x0 = x_;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute solution from x0
            obj.solve(x0, ref);
            
            % Evaluate u0
            nlp_x = obj.sol.x;
            id = obj.idx.u0;
            u = full( nlp_x(id(1):id(2)) );      
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % Delay compensation: Save current u
            if obj.expected_delay > 0
               obj.mem_u = [obj.mem_u(:,2:end), u];
            end
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargout > 1, T_opt = obj.get_T_opt(); end
            if nargout > 2, X_opt = obj.get_X_opt(); end
            if nargout > 3, U_opt = obj.get_U_opt(); end
            return
            
            % Additional evaluation
            % Complete trajectory
            % % X_opt = full(reshape(nlp_x(idx_X(1):idx_X(2)), obj.nx, obj.N));
            % % U_opt = full(reshape(nlp_x(idx_U(1):idx_U(2)), obj.nu, obj.N - 1));
            % %
            % % cost_opt = full(sol.f);
            % % constr_opt = full(sol.g);
            % %
            % % stats = obj.solver.stats;
        end
        
        function solve(obj, x0, ref)
            
            % ---- Set the initial state and reference ----
            obj.nlp_p = [x0; ref];     % Initial condition
            obj.nlp_x0(1:obj.nx) = x0; % Initial guess consistent
            
            % ---- Solve the optimization problem ----
            args = {'x0', obj.nlp_x0, ...
                'lbg', obj.nlp_lbg, ...
                'ubg', obj.nlp_ubg, ...
                'lbx', obj.nlp_lbx, ...
                'ubx', obj.nlp_ubx, ...
                'p', obj.nlp_p, ...
                %                 'lam_x0', obj.nlp_lam_x0, ...
                %                 'lam_g0', obj.nlp_lam_g0
                };
            
            obj.sol = obj.solver(args{:});
            if obj.solver.stats.success ~= true
                solve_status_str = obj.solver.stats.return_status;
                fprintf([' [' class(obj) ': ' solve_status_str '] ']);
                obj.sol.x(obj.idx.u0) = nan;
            end
            
            % Use the current solution to speed up the next optimization
            obj.nlp_x0 = obj.sol.x;
            obj.nlp_lam_x0 = obj.sol.lam_x;
            obj.nlp_lam_g0 = obj.sol.lam_g;
        end
        function T_opt = get_T_opt(obj)
            T_opt = obj.T_opt;
        end
        function X_opt = get_X_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.X;
            X_opt = full(reshape(nlp_x(id(1):id(2)), obj.nx, obj.N));
        end
        function U_opt = get_U_opt(obj)
            nlp_x = obj.sol.x;
            id = obj.idx.U;
            U_opt = full(reshape(nlp_x(id(1):id(2)), obj.nu, obj.N - 1));
        end
    end
end

