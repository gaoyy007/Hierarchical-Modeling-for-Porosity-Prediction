function createfigure(ZData1, YData1, XData1, VertexNormals1, XData2, YData2, ZData2)
%CREATEFIGURE(ZDATA1, YDATA1, XDATA1, VERTEXNORMALS1, XDATA2, YDATA2, ZDATA2)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  VERTEXNORMALS1:  surface vertexnormals
%  XDATA2:  line xdata
%  YDATA2:  line ydata
%  ZDATA2:  line zdata

%  由 MATLAB 于 20-Feb-2021 21:27:20 自动生成

% 创建 figure
figure('Tag','Print CFTOOL to Figure',...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
    'OuterPosition',[-6.33333333333333 31 1294.66666666667 697.333333333333]);

% 创建 axes
axes1 = axes('Tag','sftool surface axes');
hold(axes1,'on');

% 创建 surface
surface('ZData',ZData1,'YData',YData1,'XData',XData1,...
    'DisplayName','untitled fit 1',...
    'VertexNormals',VertexNormals1,...
    'EdgeAlpha',0.3,...
    'CData',ZData1);

% 创建 line
line(XData2,YData2,ZData2,'DisplayName','meanb vs. x, y',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',3,...
    'Marker','o',...
    'LineStyle','none');

% 创建 xlabel
xlabel('x');

% 创建 zlabel
zlabel('meanb');

% 创建 ylabel
ylabel('y');

% 取消以下行的注释以保留坐标轴的 X 范围
% xlim(axes1,[5 2095]);
% 取消以下行的注释以保留坐标轴的 Y 范围
% ylim(axes1,[137.5 412.5]);
% 取消以下行的注释以保留坐标轴的 Z 范围
% zlim(axes1,[-0.0158 0.34123]);
view(axes1,[-13.4999999999999 47.6]);
box(axes1,'on');
grid(axes1,'on');
