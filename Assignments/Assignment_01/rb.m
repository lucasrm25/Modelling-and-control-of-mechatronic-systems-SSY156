classdef rb
    %RB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        Tm = @(r,p) [r p; zeros(1,3) 1];
        rotz = @(th) [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
        rotx = @(th) [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
        roty = @(th) [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
    end
    
    methods (Static)
        function T = Tm_DH(a,alpha,d,theta)
            translation = [a*cos(theta); a*sin(theta); d];
            rotation = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha);
                        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha);
                        0           sin(alpha)             cos(alpha)            ];
            T = rb.Tm(rotation,translation);
        end
        
        function plot(V)
            quiver3( repmat(V(1,1),3,1), repmat(V(2,1),3,1), repmat(V(3,1),3,1),...
                            V(1,2:4)'-repmat(V(1,1),3,1), V(2,2:4)'-repmat(V(2,1),3,1), V(3,2:4)'-repmat(V(3,1),3,1),...
                            'LineWidth',3, 'MaxHeadSize',0.5);
            text(V(1,2),V(2,2),V(3,2),'x');
            text(V(1,3),V(2,3),V(3,3),'y');
            text(V(1,4),V(2,4),V(3,4),'z');
        end
    end
end