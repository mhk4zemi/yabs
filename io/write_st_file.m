function write_st_file(fname, results)
    disp('write the structure file - Compatible with HAWC2')
    
    data =table2array(results);

    fid = fopen(fname,'w');
    fprintf(fid,'%s\n','#1');
    fprintf(fid,'%s\t%d\n','$1',length(table2array(results)));
    format = repmat('%4.4d\t',size(table2array(results)));format(:,end) = 'r';

    for i=1:size(data,1)
        fprintf(fid, [format(i,:),'\n'],data(i,:));
    end
    fclose(fid)
end
