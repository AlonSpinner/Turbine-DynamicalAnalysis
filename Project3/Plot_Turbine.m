function [] = Plot_Turbine(v,Fig_no)
% v  - coordinate vector {q1 q2 ... qn p1 p2 ... pn}
% Fig_no - Numeric value, figure number
v = v(:);

% Parameters
L = 65;      % m
l = 10;      % m
% symbolic parameters
syms xi chi

n = 0.5*length(v);
%% Plot the basis functions
Lv = linspace(0,L,2e3);
lv = linspace(0,l,1e3);
wPf = matlabFunction(Chebypoly(xi,n-1,0,L)+eps*xi);
uPf = matlabFunction(Chebypoly(chi,n-1,0,l)+eps*chi);
dwPf = matlabFunction(diff(Chebypoly(xi,n-1,0,L),xi));

alpha = v(1:n,1).'*dwPf(L);
figure(Fig_no);clf
plot(-v(1:n,1).'*wPf(Lv),Lv,'Linewidth',3,'Color',lines(1));
    hold all
plot(-v(1:n,1).'*wPf(L)-lv*cos(alpha)+ v(n+1:end,1).'*uPf(lv)*sin(alpha),...
    L-lv*sin(alpha)-v(n+1:end,1).'*uPf(lv)*cos(alpha),'Linewidth',3,'Color',lines(1));
plot(-v(1:n,1).'*wPf(L)-l*cos(alpha)+ v(n+1:end,1).'*uPf(l)*sin(alpha),...
    L-l*sin(alpha)-v(n+1:end,1).'*uPf(l)*cos(alpha),...
    'o','MarkerSize',15,'MarkerFaceColor',lines(1),'MarkerEdgeColor',lines(1));

ylim([0 70])
axis equal

xlabel('$x$(m)','Interpreter','Latex')
ylabel('$y$(m)','Interpreter','Latex')
set(gca,'FontSize',24)
set(gcf,'color','w')
end

