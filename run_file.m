
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
% dos([hspicepath, ' -i  ./sram6Tcell.sp'])%修改路径
 system(['C:\synopsys\Hspice_Z-2007.03\BIN\hspice_mt -mt 4 -i sram6Tcell.sp -o sram6Tcell']);
% x=loadsig('./sram.tr0');%输出的bit,nbit都是10*dim的矩阵，10是由什么决定的？由网表中（tran 1ps 30ps）可知测了30组数据。有没可能是写入第3,6,9....30的数据？？最后一行是不是就是最终放电的结果？其中ndata的初始值不是0；data初始值有部分偏离1.3，也许是loadsig.c的问题
% bit=evalsig(x,'bit');%以下两行就是将x中的bit，nbit读出来
% nbit = evalsig(x,'nbit');
% 
% % - calculate diff-
% diff = bit - nbit;
% output = 1.3 - diff(end,:);%网表中设置的vdd就是1.3v。由于bit始终为1.3不变，因此输出的就是nbit的最后一行
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
while ~feof(fid) && cnt<num %feof若文件结束返回非零值，未结束返回0
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
%前面2行line代表是第2行；两段strtok代表取第2列。如果需要改变行数，就调节前面和后面的line，但总数量要保持不变。
fclose(fid);

end