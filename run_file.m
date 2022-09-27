
function  [output] = run_file(hspicepath, sample, num)
fp1 = fopen('.\sweep_data_mc','w');
fprintf(fp1,'.DATA data\n');
fprintf(fp1,'voff_pd1 vth0_pd1 ub_pd1 voff_wl1 vth0_wl1 ub_wl1\n');
for i=1:num
    fprintf(fp1, '%e\t%e\t%e\t%e\t%e\t%e\n', sample(i,1), sample(i,2), sample(i,3), sample(i,4), sample(i,5), sample(i,6));
end
fprintf(fp1,'.ENDDATA\n');
fclose(fp1);

% -simulate -
% dos([hspicepath, ' -i  ./sram6Tcell.sp'])%�޸�·��
 system(['C:\synopsys\Hspice_Z-2007.03\BIN\hspice_mt -mt 4 -i sram6Tcell.sp -o sram6Tcell']);
% x=loadsig('./sram.tr0');%�����bit,nbit����10*dim�ľ���10����ʲô�����ģ��������У�tran 1ps 30ps����֪����30�����ݡ���û������д���3,6,9....30�����ݣ������һ���ǲ��Ǿ������շŵ�Ľ��������ndata�ĳ�ʼֵ����0��data��ʼֵ�в���ƫ��1.3��Ҳ����loadsig.c������
% bit=evalsig(x,'bit');%�������о��ǽ�x�е�bit��nbit������
% nbit = evalsig(x,'nbit');
% 
% % - calculate diff-
% diff = bit - nbit;
% output = 1.3 - diff(end,:);%���������õ�vdd����1.3v������bitʼ��Ϊ1.3���䣬�������ľ���nbit�����һ��
output_origin = zeros(num,1);
fid = fopen('sram6Tcell.mt0', 'r');
%% skip the lines before alter#
while feof(fid) == 0
    line = fgetl(fid);
    if ~isempty(strfind(line,'alter#'))
        break;
    end
end

%% Get Ratio: ratio is the second string in the second row
cnt = 0;
while ~feof(fid) && cnt<num %feof���ļ��������ط���ֵ��δ��������0
    line = fgetl(fid);
    line = fgetl(fid);
    line = fgetl(fid);
%     line = fgetl(fid);
%     line = fgetl(fid);
%     line = fgetl(fid);
%     line = fgetl(fid);
    remainder = line;
    [chopped, remainder] = strtok(remainder);
if (strcmp(chopped,'failed'))
      chopped='1e-7';
end
    cnt = cnt + 1;
    output_origin(cnt,1) = str2num(chopped);
    line = fgetl(fid);
end
output=output_origin;
%ǰ��2��line�����ǵ�2�У�����strtok����ȡ��2�С������Ҫ�ı��������͵���ǰ��ͺ����line����������Ҫ���ֲ��䡣
fclose(fid);

end