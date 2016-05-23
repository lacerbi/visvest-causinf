function y = simpson1(x,dim)
%SIMPSON1 Integration with Simpson's rule along one direction.

if nargin < 1 || isempty(dim); dim = find(size(x) > 1,1); end

if dim > 2 || ~ismatrix(x)
    error('SIMPSON1 only supports 2-D matrices.')
end

w = [17 59 43 49]/48;

if dim == 1
    y = w(1)*x(1,:) + w(2)*x(2,:) + w(3)*x(3,:) + w(4)*x(4,:) + ...
        sum(x(5:end-4,:),1) + ...
        w(4)*x(end-3,:) + w(3)*x(end-2,:) + w(2)*x(end-1,:) + w(1)*x(end,:);            
else
    y = w(1)*x(:,1) + w(2)*x(:,2) + w(3)*x(:,3) + w(4)*x(:,4) + ...
        sum(x(:,5:end-4),2) + ...
        w(4)*x(:,end-3) + w(3)*x(:,end-2) + w(2)*x(:,end-1) + w(1)*x(:,end);            
end