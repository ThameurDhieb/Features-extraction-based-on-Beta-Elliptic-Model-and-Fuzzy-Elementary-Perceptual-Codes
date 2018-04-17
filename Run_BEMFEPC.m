%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%
%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  BEMFEPC  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%
%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%

Number_of_lines_simplified_visual_form = 2;
Number_of_lines_overlapped_visual_form = 4;
for Writer_Number = 1 : 1
    t = num2str(Writer_Number); 
    
    path_folder_Samples = ['..\data\',t,'\'];
    path_folder_Results = ['..\results\'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%

character_Writer_Number = num2str(Writer_Number); 
mkdir(path_folder_Results , ['writer_' , character_Writer_Number, '_AOBS']);
path_folder_Results_AOBS = [ path_folder_Results , 'writer_' , character_Writer_Number , '_AOBS', '\'];
mkdir(path_folder_Results , ['writer_' , character_Writer_Number, '_SBS']);
path_folder_Results_SBS = [ path_folder_Results , 'writer_' , character_Writer_Number, '_SBS' , '\'];
mkdir(path_folder_Results , ['writer_' , character_Writer_Number, '_AOBSFEPC']);
path_folder_Results_AOBSFEPC = [ path_folder_Results , 'writer_' , character_Writer_Number , '_AOBSFEPC', '\'];
mkdir(path_folder_Results , ['writer_' , character_Writer_Number, '_SBSFEPC']);
path_folder_Results_SBSFEPC = [ path_folder_Results , 'writer_' , character_Writer_Number , '_SBSFEPC', '\'];
delete_trace_excedent = 'oui';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

word_number = 0;
av_files = dir( fullfile(path_folder_Samples , '*.inkml') );
File_number = size(av_files,1);
matrix_param_1 = [];
matrix_param_2 = [];

for fid = 1 : File_number
    File_name = av_files(fid).name;
    path_folder_inkml = [path_folder_Samples , File_name];    
    tline = path_folder_inkml;
    [data] = lecture_online(path_folder_inkml);
    j=findstr(File_name , '.inkml');
    File_name_without_ext = File_name( 1 : j-1 );   
    word_number = word_number + 1;  
% %--------------------------------------------------------------------------
       
    [mmatrix_param_1, mmatrix_param_2] = Lecture_traitement_phrase_online_beta_elliptique_preclass_2(data, tline, word_number, Writer_Number, path_folder_Results_AOBS, path_folder_Results_AOBSFEPC, path_folder_Results_SBS, path_folder_Results_SBSFEPC, File_name_without_ext, delete_trace_excedent, Number_of_lines_simplified_visual_form, Number_of_lines_overlapped_visual_form);
    matrix_param_1 = [matrix_param_1, mmatrix_param_1];
    matrix_param_2 = [matrix_param_2, mmatrix_param_2];

    pause(0.05); 
    
end
end

