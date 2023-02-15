function save_measurement(pid,...
                          ear,...
                          mtype,...
                          data,...
                          datavars,...
                          config,...
                          configvars)

    default_save_folder = '../../DATA/'; % Must have / at the end

    time = datestr(now,'dd-mm-yyyy HH-MM-SS');
    
    filename = measnaming(pid,mtype,ear,time);

    
%     dlen = length(data);
%     clen = length(config);
    
    Config = cell2struct(config,configvars,2);
    Data = cell2struct(data,datavars,2);
    
    uisave({'Config','Data'},[default_save_folder filename])
    
%     uisave({'Config','Data'},[default_save_folder filename])

    
    % RUN THIS TO TEST:
    % A=1;B=2;C=3;D=4;E=5;pid='F04';ear='L';mtype='TE';save_measurement(pid,ear,mtype,{A,B},{'A','B'},{C,D,E},{'C','D','E'});load(uigetfile());display(Config);display(Data);
    
    

end