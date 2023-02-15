function res = measnaming(pid,mtype,ear,date)

    pidregex   = '^[FM]\d{2}$';
    mtyperegex = '^(TE|DP)$';
    earregex   = '^[LR]$';
    if regexp(pid,pidregex) ~= 1 && ...
       regexp(mtype,mtyperegex) ~= 1 && ...
       regexp(ear,earregex) ~= 1
        res = -1;
        error();
    end
    
    res = [pid ' ' mtype ' ' ear ' ' date '.mat'];
    
end