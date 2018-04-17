function [m] = lecture_online(chemin_acces)
%file='1231874526312.inkml';
%fid1=fopen(file,'rt');

%chemin_acces = 'C:\Program Files\MATLAB704\Copie de Nouveau work_13_2_Abd_AL_KARIM\seta_1\1231874526312.inkml';
%chemin_acces = 'C:\Program Files\MATLAB704\Copie de Nouveau work_13_2_Abd_AL_KARIM\seta_1\1231874543519.inkml';
%chemin_acces = 'C:\Program Files\MATLAB704\Copie de Nouveau work_13_2_Abd_AL_KARIM\seta_1\1231874560445.inkml';
%chemin_acces = 'C:\Program Files\MATLAB704\Copie de Nouveau work_13_2_Abd_AL_KARIM\seta_1\1231874595748.inkml';
%chemin_acces = 'C:\Program Files\MATLAB704\Copie de Nouveau work_13_2_Abd_AL_KARIM\seta_1\1231874575499.inkml';

fid1=fopen(chemin_acces,'rt');

m=[];
while feof(fid1) == 0

    tline = fgetl(fid1);
    if findstr(tline,'trace')> 0
        j=findstr(tline,'">');
        i=findstr(tline,'</');
        m=[m ,tline(j+2:i-1)];
        m=[m ,',0 0,','0 0,'];
        
    end
end
m=strrep(m,',',';');

m=str2num(m);

%plot(m(:,1),m(:,2),'b.');

fclose (fid1);

