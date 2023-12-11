function [ filelist ] = getfilelist( prod_type,folder_name,yyyy,startmm,startdd,endmm,enddd )
% The function could be applied to get the filelist of the .nc files for
% the user-defined period.
% User should select the folder that contains the .nc files.
% The period input must be a number (num), e.g., 2020 for year of 2020 and
% 01 for January, and 10 for October.
eval(strcat('filelist=dir(fullfile(''',folder_name,'\S5P_OFFL_L2__',...
    prod_type,'____',num2str(yyyy),num2str(startmm,'%02d'),...
    num2str(startdd,'%02d'),'*.nc''));'));
for i = startmm+1:endmm-1
    eval(strcat('filelist_',num2str(i),'=dir(fullfile(''',folder_name,...
        '\S5P_OFFL_L2__',prod_type,'____',num2str(yyyy),...
        num2str(i,'%02d'),'*.nc''));'));
    eval(strcat('b=length(filelist_',num2str(i),');'));
    eval(strcat('filelist(end+(1:b))=filelist_',num2str(i),';'));
end
switch startmm
    case {1,3,5,7,8,10,12}
        count = 31;
    case {4,6,9,11}
        count = 30;
    case 2
        count = 28;
end
diff = endmm - startmm;
switch diff
    case 0
        for j = startdd+1:enddd
            eval(strcat('filelist_',num2str(startmm),num2str(j),...
                '=dir(fullfile(''',folder_name,'\S5P_OFFL_L2__',...
                prod_type,'____',num2str(yyyy),num2str(startmm,'%02d'),...
                num2str(j,'%02d'),'*.nc''));'));
            eval(strcat('b=length(filelist_',num2str(startmm),num2str(j),');'));
            eval(strcat('filelist(end+(1:b))=filelist_',...
                num2str(startmm),num2str(j),';'));
        end
    case {1,2,3,4,5,6,7,8,9,10,11}
        for j = startdd+1:count
            eval(strcat('filelist_',num2str(startmm),num2str(j),...
                '=dir(fullfile(''',folder_name,'\S5P_OFFL_L2__',...
                prod_type,'____',num2str(yyyy),num2str(startmm,'%02d'),...
                num2str(j,'%02d'),'*.nc''));'));
            eval(strcat('b=length(filelist_',num2str(startmm),num2str(j),');'));
            eval(strcat('filelist(end+(1:b))=filelist_',...
                num2str(startmm),num2str(j),';'));
        end
        for jj = 1:enddd
            eval(strcat('filelist_',num2str(endmm),num2str(jj),...
                '=dir(fullfile(''',folder_name,'\S5P_OFFL_L2__',...
                prod_type,'____',num2str(yyyy),num2str(endmm,'%02d'),...
                num2str(jj,'%02d'),'*.nc''));'));
            eval(strcat('b=length(filelist_',num2str(endmm),num2str(jj),');'));
            eval(strcat('filelist(end+(1:b))=filelist_',...
                num2str(endmm),num2str(jj),';'));
        end
end
numfile = length(filelist);
for k = 1:numfile
    eval(strcat('[filelist(',num2str(k),...
        ').name] = [''',folder_name,'\'',filelist(',num2str(k),').name];'));
end
end