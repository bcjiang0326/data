function DataFormat(data,formatSpec,os)
% 20140714
% �����������ָ����ʽ, ��������latex �ı��
% Input:
%       data: ����
%       formatSpec: data ��һ�������ʽ string
%       os: �ļ������
N = size(data,1);
if nargin<3
    for i = 1:N
        fprintf(formatSpec,data(i,:));
    end
else
    for i = 1:N
        fprintf(os,formatSpec,data(i,:));
    end
end
end
%{
% ����
formatSpec = '%d  &  %d  &  %.6f  &  %.6f  &  %.2f  &  %.6f  &  %.6f  & %.2f  &  %.6f  &  %.6f  &  %.2f  & \\\\ \n';
data = [24,12,0.193354,0.192741  2261.36,0.193833,0.192966,2256.01,0.190243,0.188761,2443.60];
os = fopen('out.txt','w');
DataFormat(data,formatSpec,os);
DataFormat(data,formatSpec);
%}