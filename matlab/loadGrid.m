function ax = loadGrid(B, varargin)

error(nargchk(1,inf,nargin,'struct'));


[n, m] = size(B); 
if m == 6
    
    for i = 1:n
        P0 = B(i, 1:m/2);
        P1 = B(i, (m/2+1):m);
        X = [P0(1), P1(1)];
        Y = [P0(2), P1(2)];
        Z = [P0(3), P1(3)];
        
%         if (Y(1) == 1)
%             if (X(1) ~= 1 && Z(1) ~= 1)
%                 continue;
%             end
%         end
%         
%         if (X(2) == -1)
%             if (Y(2) ~= 0 && Z(1) ~= 1)
%                 continue;
%             end
%         end
%         
%         if (Z(2) == -1)
%             if (Y(2) ~= 0 && (X(1) ~= Y(1)))
%                 continue;
%             end
%         end
        
        plot3(X, Y, Z, varargin{1:end});
        hold on;
    end % for
    
else 
    
    for i = 1:n
        P0 = B(i, 1:m/2);
        P1 = B(i, (m/2+1):m);
        
        X = [P0(1), P1(1)];
        Y = [P0(2), P1(2)];
        
        plot(X, Y, varargin{1:end});
%         plot(Y, X, varargin{1:end});
        hold on; 
    end % for
    
end % if 
hold off;
ax = gca;
end % function