% Developer: Lucas Rath (https://github.com/lucasrm25)

classdef rb
    % Robotics library  
    
    properties (Constant)
        % Homogeneous transformation matrix
        Tm = @(r,p) [r p; zeros(1,3) 1];
        
        % Rotation matrices (can be used with symbolic variables)
        rotz = @(th) [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
        rotx = @(th) [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
        roty = @(th) [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
    end
    
    methods (Static)
        
%         function trajectoryPlanning(traj_table)
%         end
        
        
        % Transformation matrix that converts analytical Jacobian (ZYZ
        % Euler) to geometrical one
        % Input:    ZYZ Euler angles
        % Output:   Jacobian matrix
        %           First 3 rows = linear velocities
        %           Last  3 rows = angular velocities
        function T = T_zyz2geom(eulZYZ)
            phi   = eulZYZ(1);
            theta = eulZYZ(2);
            psi   = eulZYZ(3); %#ok<NASGU>
            Tphi = [0 -sin(phi) cos(phi)*sin(theta) ;
                    0  cos(phi) sin(phi)*sin(theta) ;
                    1  0        cos(theta)           ];
            T = [eye(3)   zeros(3);
                 zeros(3) Tphi];
        end
        
        % Calculate Geometrical Jacobian for revolute joints
        % Input:    Cell array of Homogeneous Transformation Matrices for link
        %           1 until link N
        function J_geom = J_geom_rev(T_DH_cell)
            n_links = size(T_DH_cell,1);
            
            T_0_i{1} = T_DH_cell{1};
            z{1} = [0 0 1]';                        % z0
            p{1} = [0 0 0 1]';                      % p0
            for i=2:n_links
                T_0_i{i} = T_0_i{i-1} * T_DH_cell{i};
                z{i} = T_0_i{i-1}(1:3,1:3) * z{1};    % z{2} = z1 = R_0_1 * z0
                p{i} = T_0_i{i-1}          * p{1};    % p{2} = p1 = A_0_1 * p0
            end
            pe = T_0_i{n_links} * p{1};
            
            sel = [eye(3) zeros(3,1)];              % remove last row from vector
            J_geom = sym(zeros(6,n_links-1));
            for i=1:n_links
                J_geom(1:3,i) = cross(z{i},sel*(pe-p{i}));
                J_geom(4:6,i) = z{i};
            end
            J_geom = simplify(J_geom);
        end
        
        % Homogeneous transformation matrix using Denavit-Hartenberg (DH) Convention
        % Output:   Cell array of Homogeneous Transformation Matrices for N
        %           links   % T{i} = T^{i-1}_i
        % Input:    DH_table (Nx4)  Link_i= [a_i,alpha_i,d_i,theta_i]
        function T = T_DH_table(DH_table)
            T = cell(size(DH_table,1),1);
            for i=1:size(DH_table,1)
                T{i} = rb.T_DH( DH_table(i,1),DH_table(i,2),DH_table(i,3),DH_table(i,4) );
            end
        end
        
        % Homogeneous transformation matrix using Denavit-Hartenberg Convention
        function T = T_DH(a,alpha,d,theta)
            translation = [a*cos(theta); a*sin(theta); d];
            rotation = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha);
                        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha);
                        0           sin(alpha)             cos(alpha)            ];
            T = rb.Tm(rotation,translation);
        end
        
        % Plot axis given homogeneous transformation matrix T
        function plot(T)
            quiver3( repmat(T(1,1),3,1), repmat(T(2,1),3,1), repmat(T(3,1),3,1),...
                            T(1,2:4)'-repmat(T(1,1),3,1), T(2,2:4)'-repmat(T(2,1),3,1), T(3,2:4)'-repmat(T(3,1),3,1),...
                            'LineWidth',3, 'MaxHeadSize',0.5);
            text(T(1,2),T(2,2),T(3,2),'x');
            text(T(1,3),T(2,3),T(3,3),'y');
            text(T(1,4),T(2,4),T(3,4),'z');
        end
    end
end