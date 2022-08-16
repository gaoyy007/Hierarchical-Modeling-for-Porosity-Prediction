function createfigure(ZData1, YData1, XData1, VertexNormals1, XData2, YData2, ZData2)
%CREATEFIGURE(ZDATA1, YDATA1, XDATA1, VERTEXNORMALS1, XDATA2, YDATA2, ZDATA2)
%  ZDATA1:  surface zdata
%  YDATA1:  surface ydata
%  XDATA1:  surface xdata
%  VERTEXNORMALS1:  surface vertexnormals
%  XDATA2:  line xdata
%  YDATA2:  line ydata
%  ZDATA2:  line zdata

%  �� MATLAB �� 20-Feb-2021 21:27:20 �Զ�����

% ���� figure
figure('Tag','Print CFTOOL to Figure',...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
    'OuterPosition',[-6.33333333333333 31 1294.66666666667 697.333333333333]);

% ���� axes
axes1 = axes('Tag','sftool surface axes');
hold(axes1,'on');

% ���� surface
surface('ZData',ZData1,'YData',YData1,'XData',XData1,...
    'DisplayName','untitled fit 1',...
    'VertexNormals',VertexNormals1,...
    'EdgeAlpha',0.3,...
    'CData',ZData1);

% ���� line
line(XData2,YData2,ZData2,'DisplayName','meanb vs. x, y',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',3,...
    'Marker','o',...
    'LineStyle','none');

% ���� xlabel
xlabel('x');

% ���� zlabel
zlabel('meanb');

% ���� ylabel
ylabel('y');

% ȡ�������е�ע���Ա���������� X ��Χ
% xlim(axes1,[5 2095]);
% ȡ�������е�ע���Ա���������� Y ��Χ
% ylim(axes1,[137.5 412.5]);
% ȡ�������е�ע���Ա���������� Z ��Χ
% zlim(axes1,[-0.0158 0.34123]);
view(axes1,[-13.4999999999999 47.6]);
box(axes1,'on');
grid(axes1,'on');
