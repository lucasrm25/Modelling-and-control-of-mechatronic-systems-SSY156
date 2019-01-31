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
        
        % Euler angle (ZYZ) velocities to base angular velocities
        % Analytical to geometrical Jacobian
        %   First 3 rows = linear velocities
        %   Last  3 rows = angular velocities
        function T = T_anal2geom(anal)
            phi = anal(1);
            theta = anal(2);
            % psi = anal(3);
            Tphi = [0 -sin(phi) cos(phi)*sin(theta) ;
                    0  cos(phi) sin(phi)*sin(theta) ;
                    1  0        cos(theta)           ];
            T = [eye(3)   zeros(3);
                 zeros(3) Tphi];
        end
        
        % Homogeneous transformation matrix using Denavit-Hartenberg (DH) Convention
        % Output:   cell array of T matrices given a DH table
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