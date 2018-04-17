
function [mat_CPE_p_f] = ExtractCPEF(x,grande_mmatrice_param_type1_entree_scripteurs_successifs)
xx=[];
facteur_modulo_2pi = fix( x / (2*pi) );
teta_restrainte = x - (2 * pi * facteur_modulo_2pi);
zise_vect_teta = size(teta_restrainte , 1);
for k_Th_Dh = 1 : zise_vect_teta
    if ( teta_restrainte(k_Th_Dh , 1) > pi )
        teta_restrainte(k_Th_Dh , 1) = teta_restrainte(k_Th_Dh , 1) - (2*pi);
     elseif ( teta_restrainte(k_Th_Dh , 1) < -pi )
       teta_restrainte(k_Th_Dh , 1) = teta_restrainte(k_Th_Dh , 1) + (2*pi);
    end
end
x=teta_restrainte;
[M5,N5]=size(x);
if x(M5,1)==0
    x(M5,1) = 0.01;
end



% Création des contrôleurs flou Ci de type Sugeno

% FIS=NEWFIS(FISNAME, FISTYPE) creates a FIS structure for a Mamdani or 
% Sugeno-style system with the name FISNAME.
 
C1=newfis('C1', 'sugeno');
C2=newfis('C2', 'sugeno');
C3=newfis('C3', 'sugeno');
C4=newfis('C4', 'sugeno');
C5=newfis('C5', 'sugeno');
C6=newfis('C6', 'sugeno');
C7=newfis('C7', 'sugeno');
C8=newfis('C8', 'sugeno');

% Construction des sous ensembles flous associés "input"  
% %%% valeurs de theta de l'intervalle -pi : 0 : pi

%%%% définition des différentes bornes des parties de thétas
b1_p1= -pi;         %%% borne  1 de la partie N°1
b2_p1 = -7*pi/8;    %%% borne  2 de la partie N°1 

b1_p2 = -7*pi/8;   %%% borne  1 de la partie N°2
b2_p2 = -5*pi/8;   %%% borne  2 de la partie N°2 

b1_p3 = -5*pi/8;   %%% borne  1 de la partie N°3
b2_p3 = -3*pi/8;   %%% borne  2 de la partie N°3 

b1_p4 = -3*pi/8;   %%% borne  1 de la partie N°4
b2_p4 = -pi/8;     %%% borne  2 de la partie N°4

b1_p5 = -pi/8;   %%% borne  1 de la partie N°5
b2_p5 = pi/8;    %%% borne  2 de la partie N°5 

b1_p6 = pi/8;   %%% borne  1 de la partie N°6
b2_p6 = 3*pi/8; %%% borne  2 de la partie N°6

b1_p7 = 3*pi/8;  %%% borne  1 de la partie N°7
b2_p7 = 5*pi/8;  %%% borne  2 de la partie N°7 

b1_p8 = 5*pi/8;  %%% borne  1 de la partie N°8
b2_p8 = 7*pi/8;  %%% borne  2 de la partie N°8 

b1_p9 = 7*pi/8;  %%% borne  1 de la partie N°9
b2_p9 = pi;      %%% borne  2 de la partie N°9 
%%%%%

%%%%%%
const=pi/16;%%%% valeur de chauvechement en radian ( 11.25°)
%const=pi/12;  %%% 15 degres
% const=pi/8;  %%%% 22.5 °
% const=pi/32;
const_degre=const*180/pi; %%%% valeur de chauvechement en degree

% Création des variables d entrées et de sortie du contrôleur flou
C1=addvar(C1,'input','theta',[-pi pi]);
C1=addvar(C1,'output','EPC',[0 2]);
C2=addvar(C2,'input','theta',[-pi pi]);
C2=addvar(C2,'output','EPC',[0 2]);
C3=addvar(C3,'input','theta',[-pi pi]);
C3=addvar(C3,'output','EPC',[0 2]);
C4=addvar(C4,'input','theta',[-pi pi]);
C4=addvar(C4,'output','EPC',[0 2]);
C5=addvar(C5,'input','theta',[-pi pi]);
C5=addvar(C5,'output','EPC',[0 2]);
C6=addvar(C6,'input','theta',[-pi pi]);
C6=addvar(C6,'output','EPC',[0 2]);
C7=addvar(C7,'input','theta',[-pi pi]);
C7=addvar(C7,'output','EPC',[0 2]);
C8=addvar(C8,'input','theta',[-pi pi]);
C8=addvar(C8,'output','EPC',[0 2]);
% intervalle_a_c_tri=b2_p2+const-((b1_p2-const)+(b2_p2+const))/2;%%%%%valeur de l'intervalle b-c de la fonction tringulaire en radian
% intervalle_a_c_tri_degree=intervalle_a_c_tri*180/pi; %%%%%valeur de l'intervalle b-c de la fonction triangulaire en degree

C1=addmf(C1,'input',1,'EPN','trimf',[b1_p1-1  b1_p1 (b1_p2)+const]);  
C1=addmf(C1,'input',1,'NP','trimf',[b1_p2-const ((b1_p2-const)+(b2_p2+const))/2  b2_p2+const ]);

C2=addmf(C2,'input',1,'SN','trimf',[b1_p2-const ((b1_p2-const)+(b2_p2+const))/2  b2_p2+const ]);
C2=addmf(C2,'input',1,'MN','trimf',[b1_p3-const ((b1_p3-const)+(b2_p3+const))/2  b2_p3+const]);

C3=addmf(C3,'input',1,'NM','trimf',[b1_p3-const ((b1_p3-const)+(b2_p3+const))/2  b2_p3+const]);
C3=addmf(C3,'input',1,'NG','trimf',[b1_p4-const ((b1_p4-const)+(b2_p4+const))/2  b2_p4+const]);

C4=addmf(C4,'input',1,'NG','trimf',[b1_p4-const ((b1_p4-const)+(b2_p4+const))/2  b2_p4+const]);
C4=addmf(C4,'input',1,'EZ','trimf',[b1_p5-const ((b1_p5-const)+(b2_p5+const))/2  b2_p5+const]);

C5=addmf(C5,'input',1,'EZ','trimf',[b1_p5-const ((b1_p5-const)+(b2_p5+const))/2  b2_p5+const]);
C5=addmf(C5,'input',1,'PP','trimf',[b1_p6-const ((b1_p6-const)+(b2_p6+const))/2  b2_p6+const]);

C6=addmf(C6,'input',1,'PP','trimf',[b1_p6-const ((b1_p6-const)+(b2_p6+const))/2  b2_p6+const]);
C6=addmf(C6,'input',1,'PM','trimf',[b1_p7-const ((b1_p7-const)+(b2_p7+const))/2  b2_p7+const]);

C7=addmf(C7,'input',1,'PM','trimf',[b1_p7-const ((b1_p7-const)+(b2_p7+const))/2  b2_p7+const]);
C7=addmf(C7,'input',1,'PG','trimf',[b1_p8-const ((b1_p8-const)+(b2_p8+const))/2  b2_p8+const]);

C8=addmf(C8,'input',1,'PG','trimf',[b1_p8-const ((b1_p8-const)+(b2_p8+const))/2  b2_p8+const]);
C8=addmf(C8,'input',1,'EPP','trimf',[b1_p9-const b2_p9  b2_p9+1]);

% % Représentation graphique des sous ensembles flous  "input"
% figure(1)
% plotmf(C1,'input',1)

% % Représentation graphique des sous ensembles flous  "output"

C1=addmf(C1,'output',1,'Vallee','constant',1);
C1=addmf(C1,'output',1,'H-O-G','constant',2);
%
C2=addmf(C2,'output',1,'Left oblique shaft','constant',1);
C2=addmf(C2,'output',1,'Shaft','constant',2);
%
C3=addmf(C3,'output',1,'Hampe','constant',1);
C3=addmf(C3,'output',1,'H-O-D','constant',2);
%
C4=addmf(C4,'output',1,'H-O-D','constant',1);
C4=addmf(C4,'output',1,'Vallee','constant',2);
%
C5=addmf(C5,'output',1,'Vallee','constant',1);
C5=addmf(C5,'output',1,'H-O-G','constant',2);
%
C6=addmf(C6,'output',1,'H-O-G','constant',1);
C6=addmf(C6,'output',1,'Hampe','constant',2);
%
C7=addmf(C7,'output',1,'Hampe','constant',1);
C7=addmf(C7,'output',1,'H-O-D','constant',2);
%
C8=addmf(C8,'output',1,'H-O-D','constant',1);
C8=addmf(C8,'output',1,'Vallee','constant',2);
%%%%%
ruleList1=[ ...
1 1 1 2
2 2 1 2
];
C1=addrule(C1,ruleList1);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C1)
% fuzzy(C1)
% pause
%%
ruleList2=[ ...
1 1 1 2
2 2 1 2
];
C2=addrule(C2,ruleList2);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C2)
% fuzzy(C2)
%%
ruleList3=[ ...
1 1 1 2
2 2 1 2
];
C3=addrule(C3,ruleList3);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C3)
% fuzzy(C3)
%%
ruleList4=[ ...
1 1 1 2
2 2 1 2
];
C4=addrule(C4,ruleList4);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C4)
% fuzzy(C4)
%%
ruleList5=[ ...
1 1 1 2
2 2 1 2
];
C5=addrule(C5,ruleList5)
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C5)
% fuzzy(C5)
%%
ruleList6=[ ...
1 1 1 2
2 2 1 2
];
C6=addrule(C6,ruleList6);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C6)
% fuzzy(C6)
%%
ruleList7=[ ...
1 1 1 2
2 2 1 2
];
C7=addrule(C7,ruleList7);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C7)
% fuzzy(C7)
%%
ruleList8=[ ...
1 1 1 2
2 2 1 2
];
C8=addrule(C8,ruleList8);
% x=pi/16;
% x_degre=(x*180)/pi
% out=evalfis(x,C8)
% fuzzy(C8)

%%% calcul du pourcentage d'appartenance  de chaque cpe
%%% index n°1: CPE1 : poucentage d'appartenance au CPE1(VALLEE)
%%% index n°2: CPE2 : poucentage d'appartenance au CPE2(H-O-G)
%%% index n°3: CPE3 : poucentage d'appartenance au CPE3(hampe)
%%% index n°4: CPE4 : poucentage d'appartenance au CPE4(H-O-D)
% 
% load result_flou % du fichier exe_tablet_ext_rad.m line 303
% load('DebutHampeDebut_ICEP.mat');
mat_teta_f_pi=x;
% mat_teta_f_pi
[M,N]=size(mat_teta_f_pi);
% mat_CPE_p=[];
mat_out=[];
%%% evaluation du système pour le controleur flou C1 
%%% les CPE vallée et H-O-G
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p1-1)&(mat_teta_f_pi(i,j)<=b1_p3-const ))

       
out=evalfis([mat_teta_f_pi(i,j)],C1); %%% evalution du système    
          mat_out(i,j)=out;
        %%% CPE1 vallee
        if (out==1)
            mat_CPE1(i,j)=100;
            mat_CPE2(i,j)=0;            
            mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE1(i,j)=100-((abs(1-out))*100);
            mat_CPE2(i,j)=abs(1-out)*100;
            mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out>1.5)
            mat_CPE2(i,j)=abs(1-out)*100; 
            mat_CPE1(i,j)=100-((abs(1-out))*100);
            mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        %%% CPE2 H-O-G
        if (out==2)
            mat_CPE2(i,j)=100; 
            mat_CPE1(i,j)=0;
            mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
    end
end
end
%%% evaluation du système pour le controleur flou C2 
%%% les CPE H-O-G et Hampe
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p3-const)&(mat_teta_f_pi(i,j)<=b1_p4-const ))
        out=evalfis([mat_teta_f_pi(i,j)],C2); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE2 H-O-G
        if (out==1)
            mat_CPE2(i,j)=100;
            mat_CPE3(i,j)=0;  
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE2(i,j)=100-((abs(1-out))*100);
            mat_CPE3(i,j)=abs(1-out)*100;
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out>1.5)
            mat_CPE3(i,j)=abs(1-out)*100; 
            mat_CPE2(i,j)=100-((abs(1-out))*100);
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        %%% CPE Hampe
        if (out==2)
            mat_CPE3(i,j)=100; 
            mat_CPE2(i,j)=0;
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
    end
end
end
%%% evaluation du système pour le controleur flou C3 
%%% les CPE Hampe et H-O-D
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p4-const)&(mat_teta_f_pi(i,j)<=b1_p5-const ))
        out=evalfis([mat_teta_f_pi(i,j)],C3); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE3 Hampe
        if (out==1)
            mat_CPE3(i,j)=100;
            mat_CPE4(i,j)=0; 
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE3(i,j)=100-((abs(1-out))*100);
            mat_CPE4(i,j)=abs(1-out)*100;
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
        if (out>1.5)
            mat_CPE4(i,j)=abs(1-out)*100; 
            mat_CPE3(i,j)=100-((abs(1-out))*100);
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
        %%% CPE H-O-D
        if (out==2)
            mat_CPE4(i,j)=100; 
            mat_CPE3(i,j)=0;
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
    end
end
end
%%% evaluation du système pour le controleur flou C4
%%% les CPE H-O-D et vallée
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p5-const)&(mat_teta_f_pi(i,j)<=b1_p6-const))
        out=evalfis([mat_teta_f_pi(i,j)],C4); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE4 H-O-D
        if (out==1)
            mat_CPE4(i,j)=100;
            mat_CPE1(i,j)=0;   
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE4(i,j)=100-((abs(1-out))*100);
            mat_CPE1(i,j)=abs(1-out)*100;
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
        if (out>1.5)
            mat_CPE1(i,j)=abs(1-out)*100; 
            mat_CPE4(i,j)=100-((abs(1-out))*100);
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
        %%% CPE H-O-D
        if (out==2)
            mat_CPE1(i,j)=100; 
            mat_CPE4(i,j)=0;
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
    end
end
end
%%% evaluation du système pour le controleur flou C5
%%% les CPE vallée et  H-O-G
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p6-const)&(mat_teta_f_pi(i,j)<=b1_p7-const))
        out=evalfis([mat_teta_f_pi(i,j)],C5); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE1 vallée
        if (out==1)
            mat_CPE1(i,j)=100;
            mat_CPE2(i,j)=0;  
            mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE1(i,j)=100-((abs(1-out))*100);
            mat_CPE2(i,j)=abs(1-out)*100;
            mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out>1.5)
            mat_CPE2(i,j)=abs(1-out)*100; 
            mat_CPE1(i,j)=100-((abs(1-out))*100);
             mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        %%% CPE H-O-G
        if (out==2)
            mat_CPE2(i,j)=100; 
            mat_CPE1(i,j)=0;
             mat_CPE3(i,j)=0;
            mat_CPE4(i,j)=0;
        end
    end
end
end
%%% evaluation du système pour le controleur flou C6
%%% les CPE H-O-G et hampe
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p7-const)&(mat_teta_f_pi(i,j)<=b1_p8-const ))
        out=evalfis([mat_teta_f_pi(i,j)],C6); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE2 H-O-G
        if (out==1)
            mat_CPE2(i,j)=100;
            mat_CPE3(i,j)=0; 
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE2(i,j)=100-((abs(1-out))*100);
            mat_CPE3(i,j)=abs(1-out)*100;
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        if (out>1.5)
            mat_CPE3(i,j)=abs(1-out)*100; 
            mat_CPE2(i,j)=100-((abs(1-out))*100);
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
        %%% CPE HAMPE
        if (out==2)
            mat_CPE3(i,j)=100; 
            mat_CPE2(i,j)=0;
            mat_CPE1(i,j)=0;
            mat_CPE4(i,j)=0;
        end
    end
end
end
%%% evaluation du système pour le controleur flou C7
%%% les CPE hampe et H-O-D
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p8-const)&(mat_teta_f_pi(i,j)<=b1_p9-const ))
        out=evalfis([mat_teta_f_pi(i,j)],C7); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE3 hampe
        if (out==1)
            mat_CPE3(i,j)=100;
            mat_CPE4(i,j)=0;   
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE3(i,j)=100-((abs(1-out))*100);
            mat_CPE4(i,j)=abs(1-out)*100;
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
        if (out>1.5)
            mat_CPE4(i,j)=abs(1-out)*100; 
            mat_CPE3(i,j)=100-((abs(1-out))*100);
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
        %%% CPE H-O-D
        if (out==2)
            mat_CPE4(i,j)=100; 
            mat_CPE3(i,j)=0;
            mat_CPE1(i,j)=0;
            mat_CPE2(i,j)=0;
        end
    end
end
end
%% evaluation du système pour le controleur flou C8
%%% les CPE H-O-D  et vallée
for i=1:M
     for j=1:N
        if (mat_teta_f_pi(i,j)~=0 )& ((mat_teta_f_pi(i,j)>=b1_p9-const)&(mat_teta_f_pi(i,j)<=b2_p9 ))
        out=evalfis([mat_teta_f_pi(i,j)],C8); %%% evalution du système
        mat_out(i,j)=out;
        %%% CPE4 H-O-D
        if (out==1)
            mat_CPE4(i,j)=100;
            mat_CPE1(i,j)=0;   
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
        if (out<=1.5)
            mat_CPE4(i,j)=100-((abs(1-out))*100);
            mat_CPE1(i,j)=abs(1-out)*100;
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
        if (out>1.5)
            mat_CPE1(i,j)=abs(1-out)*100; 
            mat_CPE4(i,j)=100-((abs(1-out))*100);
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
        %%% CPE vallée
        if (out==2)
            mat_CPE1(i,j)=100; 
            mat_CPE4(i,j)=0;
            mat_CPE2(i,j)=0;
            mat_CPE3(i,j)=0;
        end
    end
end
end

mat_CPE_p=[mat_CPE1 mat_CPE2 mat_CPE3 mat_CPE4];
mat_CPE_p_f=[];

for i = 1 : N
 mat_CPE_p_f = [mat_CPE_p_f mat_CPE1(:,i) mat_CPE2(:,i) mat_CPE3(:,i) mat_CPE4(:,i)];
end
% mat_CPE_p_f
mat_CPE_p_f=mat_CPE_p_f';


ALLALL=[];
xx=grande_mmatrice_param_type1_entree_scripteurs_successifs;
[n m] =size(xx);
kk=0;
for ii=1:m
% xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([1:7],:);
xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([1:8],:);
kk=kk+1;
All=[xx(:,ii);mat_CPE_p_f(:,kk)];
% xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([9:17],:);
xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([9:18],:);
kk=kk+1;
All2=[xx(:,ii);mat_CPE_p_f(:,kk)];
% xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([19:27],:);
xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([19:28],:);
kk=kk+1;
All3=[xx(:,ii);mat_CPE_p_f(:,kk)];
% xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([29:37],:);
xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([29:38],:);
kk=kk+1;
All4=[xx(:,ii);mat_CPE_p_f(:,kk)];
% xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([39:40],:);
xx=grande_mmatrice_param_type1_entree_scripteurs_successifs([39:40],:);
All=[All;All2];
All=[All;All3];
All=[All;All4];
All=[All;xx(:,ii)];
ALLALL=[ALLALL, All];
end
mat_CPE_p_f= ALLALL;

end

