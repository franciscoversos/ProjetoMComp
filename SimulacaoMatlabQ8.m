%%-------------------------------------------------------------------------
% Simulação do Problema de Escoamento, com elementos quadrangulares
% quadráticos,Q8.
%%-------------------------------------------------------------------------
%Importar os ficheiros do Nx, que contém os nós e os elementos, e outro com
%os resultados, para comparar

format long
clear all
fileID0 = fopen('geral_quad8_nodes.txt', 'r');
fileID1 = fopen('geral_quad8_elements.txt', 'r');
Q8NxPotencial= xlsread("Q8_malhabase_potencial.xlsx");        %Ficheiro com os dados do potencial do NX
Q8NxVelocidade= xlsread("Q8_malhabase_velocidade.xlsx");      %Ficheiro com os dados da velocidade do Nx
formatSpec = '%c';'%d' ;
%%-------------------------------------------------------------------------
%Retirar as coordenadas dos nós, armazenadas no vetor coord
nodes = splitlines(fscanf(fileID0, formatSpec));
auxN =[];
for i =11:6:length(nodes)
    auxN = [auxN, nodes(i)];
end
auxN = auxN';
col = split(auxN);
El1=[];
El2 =[];
coord=[];
for i =1:length(col)
    El1 = [El1; str2double(col(i,5))];
    El2 = [El2; str2double(col(i,6))];
    coord = [coord; i, El1(i), El2(i)];
end
%%-------------------------------------------------------------------------
%Criação da matriz de conectividades em que para cada elemento(Coluna1)
%correspondem-lhe os respetivos nós, armazenados no vetor Elements

Big = splitlines(fscanf(fileID1, formatSpec));
auxE =[];
dauxE=[];
for i = 17:length(Big)-4
    auxE = [auxE,Big(i)];
end
auxE = auxE';
colE = split(auxE);
El1 =[];
El2=[];
El3 =[];
El4=[];
El5 =[];
El6=[];
El7 =[];
El8=[];
Elements=[];
for i =1:length(colE)
    El1 = [El1; str2double(colE(i,11))];
    El2 = [El2; str2double(colE(i,12))];
    El3 = [El3;str2double(colE(i,13))];
    El4 = [El4; str2double(colE(i,14))];
    El5 = [El5; str2double(colE(i,15))];
    El6 = [El6; str2double(colE(i,16))];
    El7 = [El7;str2double(colE(i,17))];
    El8 = [El8; str2double(colE(i,18))];
    Elements = [Elements; i, El1(i), El2(i), El3(i),El4(i), El5(i),El6(i),El7(i),El8(i)];
end
Kr =sparse(length(coord), length(coord));
fr=sparse(length(coord),1);
%%--------------------------------------------------------------------------
%Cálculo da Matriz de Rigidez (Kr) e do vetor de forças (fr)
for i =1:length(Elements)
    no1= Elements(i,2);
    no2= Elements(i,3);
    no3= Elements(i,4);
    no4= Elements(i,5);
    no5= Elements(i,6);
    no6= Elements(i,7);
    no7= Elements(i,8);
    no8= Elements(i,9);
    
  edofs = [no1 no2 no3 no4 no5 no6 no7 no8];
  XN = [coord(no1,2), coord(no1,3);coord(no2,2), coord(no2,3);coord(no3,2), coord(no3,3);coord(no4,2), coord(no4,3); ...
      coord(no5,2), coord(no5,3);coord(no6,2), coord(no6,3);coord(no7,2), coord(no7,3);coord(no8,2), coord(no8,3)];
  f1 =0;
  [Ke fe B psi] = Elem_Quad8(XN, f1);
  Kr(edofs,edofs)=Kr(edofs,edofs)+Ke;
  fr(edofs,1)=fr(edofs,1)+fe;
end
%%-------------------------------------------------------------------------
%Condições de Fronteira: 

%Condições de Dirichlet
k = boundary(coord(:,2), coord(:,3), 0.4);
fronteira =[coord(k,1), coord(k,2), coord(k,3)];
plot(fronteira(:,2), fronteira(:,3))
coordfront= [];
for i =1:length(fronteira)
    if fronteira(i,2) == 2000
        coordfront = [coordfront;fronteira(i,1), 0];
    end 
end
boom = 10^20;
for i = 1:length(coordfront)
    Kr(coordfront(i,1),coordfront(i,1))= boom;
    fr(coordfront(i,1))= boom*coordfront(i,2);
end

%Criação dos vetores, fronteira1 que é a parede da esquerda, fronteira2 que
%é a parede de cima, fronteira3 que é a parede da direita e fronteira4 que
%é a parede de baixo
v = 2500;
fronteira1=[];
fronteira2=[];
fronteira3=[];
fronteira4=[];
n =1;
for i = 1:length(fronteira)-1
    if fronteira(i,2)==0
        fronteira1=[fronteira1;n,fronteira(i,1), fronteira(i,2), fronteira(i,3)];
        n =n+1;
    end
end
for l =1:length(fronteira)
    if fronteira(l,3)==500 
        fronteira2=[fronteira2; l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
    if fronteira(l,2)>=1000 && fronteira(l,2)<=1600 &&fronteira(l,3)>=200 && fronteira(l,3)<500
        fronteira2=[fronteira2;l, fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
end

for l=1:length(fronteira)
        if fronteira(l,2) >=1750 && fronteira(l,2)<=2000 &&fronteira(l,3)<=125
            fronteira3=[fronteira3;l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
        end
end
for l =1:length(fronteira)
    if fronteira(l,2)==2000
        fronteira4 = [fronteira4;l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
end
fronteira1 = sortrows(fronteira1,4);
fronteira2 = sortrows(fronteira2,3);
fronteira3 = sortrows(fronteira3,3);
fronteira4 = sortrows(fronteira3,3);
N = length(fronteira1);


impares1 =[];
pares1 =[];
impares2 =[];
pares2 =[];
impares3 =[];
pares3 =[];
impares4 =[];
pares4 =[];
for l = 1:N
    fronteira1(l,1)=l;
    impares1 =[impares1; 2*l-1;l];
    pares1 =[pares1; l;2*l];
end
for l = 1:length(fronteira2)
    fronteira2(l,1)=l;
    impares2 =[impares2; 2*l-1;l];
    pares2 =[pares2;1; 2*l];
end
for l = 1:length(fronteira3)
    fronteira3(l,1)=l;
    impares3 =[impares3; 2*l-1;l];
    pares3 =[pares3;1; 2*l];
end

for l = 1:length(fronteira4)
    fronteira4(l,1)=l;
    impares4 =[impares4; 2*l-1;l];
    pares4 =[pares4;1; 2*l];
end
%Condições de Neumann
%De acordo com as forças nodais, as forças distribuem-se da seguinte
%maneira
nEl=0;
n=0;
k=0;
for l =1:1:N
    if  fronteira1(l,1)== impares1(l) && (fronteira1(l,1)~=1 && fronteira1(l,1) ~=N)
        fr(fronteira1(l,2)) = (1/3)*v;
        k=k+1;
    end
    if fronteira1(l,1)== 1||fronteira1(l,1)==N
        fr(fronteira1(l,2)) = (1/6)*v;
    end
    if fronteira1(l,1)== pares1(l) && (fronteira1(l,1)~=1 &&fronteira1(l,1) ~=N)
        fr(fronteira1(l,2)) = (4/6)*v;
        n=n+1;
    end
nEl=(n+k+1)/2;
end

fr(fronteira(:,1)) =fr(fronteira(:,1))*(500/nEl);
%%-------------------------------------------------------------------------
%Cálculo do potencial por Backslash
u = Kr\fr;

%%------------------------------------------------------------------------
%Exportando o ficheiro com os dados do potencial do Nx, calcula-se o erro entre as duas

EPotencial=[];
for i =1:length(coord)
    EPotencial=[EPotencial; u(i)-Q8NxPotencial(i,8)*10^6];
end

%%-------------------------------------------------------------------------
%Cálculo dos gradientes nos pontos de interesse(SfluxuQuiver) e nos
%centróides(Sfluxu) e da pressão

L=0;
Sfluxu=[];
Xpint=[];
SfluxuQuiver=[];
for i=1:length(Elements);
no1= Elements(i,2);
    no2= Elements(i,3);
    no3= Elements(i,4);
    no4= Elements(i,5);
    no5= Elements(i,6);
    no6= Elements(i,7);
    no7= Elements(i,8);
    no8= Elements(i,9);

  edofs = [no1 no2 no3 no4 no5 no6 no7 no8];
  XN = [coord(no1,2), coord(no1,3);coord(no2,2), coord(no2,3);coord(no3,2), coord(no3,3);coord(no4,2), coord(no4,3); ...
      coord(no5,2), coord(no5,3);coord(no6,2), coord(no6,3);coord(no7,2), coord(no7,3);coord(no8,2), coord(no8,3)];
  f1 =0;
%   O centroide esta na origem
csi=0;
eta=0;
nip = 4 ;    %  pts de integracao reduzida 
[xp wp]=Genip2DQ (nip);

%   percorrer os pontos de integracao
for ip=1:nip
L=L+1 ;
csi = xp(ip,1);
eta = xp(ip,2);
%   para cada ponto calcular
%----------------------------------------------------------------
[B psi Detj]=Shape_N_Der8 (XN,csi,eta);
%----------------------------------------------------------------
uint = psi'*u(edofs);
xpint = XN'*psi;
gradu = B'*u(edofs);
fluxu = -gradu;
Xpint=[Xpint, xpint];
SfluxuQuiver = [SfluxuQuiver,fluxu];%Fluxo nos pontos de interesse, permite um mapa melhor

end

Sfluxu = [Sfluxu,fluxu];%Fluxo nos centroides
pressao(i)=(0.5).*(2500^2 -(Sfluxu(1,i).^2)+ (Sfluxu(2,i).^2));%Pressão nos centróides
end

%%-------------------------------------------------------------------------
%Cálculo das velocidades nos pontos de interesse

fno1 = 1;fno2=2;fno3 =3;fno4=4;
fno5=5;fno6=6;fno7=7;fno8=8;
fno=[fno2 fno3 fno4 fno5 fno6 fno7 fno8 fno1];
fno=fno';
W=0;
K=[];
El=[];
U=[];
%Cálculo das velocidades nos nós, fazendo a média com os valores dos
%elementos

for l =1:length(coord)
W=0;
ufno=0;
vfno=0;
    for i = 1:length(Elements)
        for k=2:9
            if coord(l) == Elements(i,k)
            W=W+1;
            ufno = ufno + Sfluxu(1,Elements(i,1));
            vfno = vfno +Sfluxu(2,Elements(i,1)) ;
            end
        end
    end
U(l,:)=[sqrt((ufno/W)^2 + (vfno/W)^2)];
end
U(fno);
u(fno);
EVelocidade=[];

%Cálculo do Erro da Velocidade entre o Matlab e o Nx
for i =1:length(coord)
    EVelocidade=[EVelocidade; U(i)-Q8NxVelocidade(i,11)*10^9];
end
%%-------------------------------------------------------------------------
%Forca resultante na fronteira 1

fronteira1d=[];
for i =1:1:length(fronteira1)
    if fronteira1(i,1) == impares1(i)
        fronteira1d =[fronteira1d; fronteira1(i,:)];
    end
end
D1=[];
f1x=[];
f1y=[];
p1=0;
p2=0;
d=0;

for l =1:length(fronteira1d)-1
    p1 =[fronteira1d(l,3); fronteira1d(l,4)];
    p2 =[fronteira1d(l+1,3);fronteira1d(l+1,4)];
    dh=p2-p1;
    d =norm(p2-p1);
    dv=[dh(2)/d;-dh(1)/d];      %Vetor da normal exterior
    D1=[D1,d];
    fx=0;
    for i=1:length(Elements)
         for k =2:5
             for y=2:5
                if fronteira1d(l,2)== Elements(i,k)&&fronteira1d(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy=pressao(Elements(i,1))*dv(2);
                end
            end
         end
    end
f1x=[f1x;fx*d];
f1y=[f1y;fy*d];
end
F1x=sum(f1x);
F1y=sum(f1y);
%%-------------------------------------------------------------------------
%Forca resultante na fronteira 2

fronteira2d=[];
for i =1:1:length(fronteira2)
    if fronteira2(i,1) == impares2(i)
        fronteira2d =[fronteira2d; fronteira2(i,:)];
    end
end

D2=[];
f2x=[];
f2y=[];
p1=0;
p2=0;
d=0;

for l =1:length(fronteira2d)-1
    p1 =[fronteira2d(l,3); fronteira2d(l,4)];
    p2 =[fronteira2d(l+1,3);fronteira2d(l+1,4)];
    dh=p2-p1;
    d =norm(p2-p1);
    dv=[dh(2)/d;-dh(1)/d];      %Vetor da normal exterior
    D2=[D2,d];
    fx=0;
    for i=1:length(Elements)
         for k =2:5
             for y=2:5
                if fronteira2d(l,2)== Elements(i,k)&&fronteira2d(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy=pressao(Elements(i,1))*dv(2);
                end
            end
         end
    end
f2x=[f2x;fx*d];
f2y=[f2y;fy*d];
end
F2x=sum(f2x);
F2y=sum(f2y);
%%-------------------------------------------------------------------------
%Forca resultante na fronteira 2

fronteira3d=[];
for i =1:1:length(fronteira3)
    if fronteira3(i,1) == impares3(i)
        fronteira3d =[fronteira3d; fronteira3(i,:)];
    end
end
D3=[];
f3x=[];
f3y=[];
p1=0;
p2=0;
d=0;

for l =1:length(fronteira3d)-1
    p1 =[fronteira3d(l,3); fronteira3d(l,4)];
    p2 =[fronteira3d(l+1,3);fronteira3d(l+1,4)];
    dh=p2-p1;
    d =norm(p2-p1);
    dv=[dh(2)/d;-dh(1)/d];      %Vetor da normal exterior
    D3=[D3,d];
    fx=0;
    for i=1:length(Elements)
         for k =2:5
             for y=2:5
                if fronteira3d(l,2)== Elements(i,k)&&fronteira3d(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy=pressao(Elements(i,1))*dv(2);
                end
            end
         end
    end
f3x=[f3x;fx*d];
f3y=[f3y;fy*d];
end
F3x=sum(f3x);
F3y=sum(f3y);
%%-------------------------------------------------------------------------
%Força resultante na fronteira4

%Como um elemento é delimitado pelos dois vértices, excluem-se os nós
%interiores(pares)
fronteira4d=[];
for i =1:1:length(fronteira4)
    if fronteira4(i,1) == impares4(i)
        fronteira4d =[fronteira4d; fronteira4(i,:)];
    end
end
D4=[];
f4x=[];
f4y=[];
p1=0;
p2=0;
d=0;

for l =1:length(fronteira4d)-1
    p1 =[fronteira4d(l,3); fronteira4d(l,4)];
    p2 =[fronteira4d(l+1,3);fronteira4d(l+1,4)];
    dh=p2-p1;
    d =norm(p2-p1);
    dv=[dh(2)/d;-dh(1)/d];      %Vetor da normal exterior
    D4=[D4,d];
    fx=0;
    for i=1:length(Elements)
         for k =2:5
             for y=2:5
                if fronteira4d(l,2)== Elements(i,k)&&fronteira4d(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy=pressao(Elements(i,1))*dv(2);
                end
            end
         end
    end
f4x=[f4x;fx*d];
f4y=[f4y;fy*d];
end
F4x=sum(f4x);
F4y=sum(f4y);

%Força resultante na peça total

Fx =(F1x +F2x +F3x+F4x)*10^-6;
Fy=(F1y + F2y +F3y + F4y)*10^-6;
%%-------------------------------------------------------------------------
%Menu que permite ao utilizador escolher que gráfico ver
while true
    % Menu 
    msg = 'Escolha o gráfico:';
    opts = {'Campo de velocidades', 'Potencial', 'Pressão','Erro Potencial','Erro Velocidade', 'Sair'};
    choice = menu(msg, opts);

    
    switch choice
        case 1 % Campo de Velocidades
            for i=1:length(Xpint)    
                plot(Xpint(1,i),Xpint(2,i),'bx')
                quiver(Xpint(1,i),Xpint(2,i),SfluxuQuiver(1,i),SfluxuQuiver(2,i), 1/20);hold on
                plot(coord(1),coord(2),'ro');
            end
            title('Campo de Velocidades');
            plot(fronteira(:,2),fronteira(:,3))
            hold off

        case 2 % Potencial
            figure;
             for i=1:length(Elements)
                no1= Elements(i,2);
                no2= Elements(i,3);
                no3= Elements(i,4);
                no4= Elements(i,5);
                no5= Elements(i,6);
                no6= Elements(i,7);
                no7= Elements(i,8);
                no8= Elements(i,9);
                edofs=[no1 no5 no2 no6 no3 no7 no4 no8]; % ordem para desenhar    
                
                fill (coord(edofs,2),coord(edofs,3),u(edofs));hold on
                plot(coord(edofs,2),coord(edofs,3),'r');hold on
                colorbar
             end
            title('Campo Potencial');
            plot(coord(1),coord(2),'ro');
            hold off
        case 3 % Pressão
            figure;
            for i=1:length(Elements)
                no1= Elements(i,2);
                no2= Elements(i,3);
                no3= Elements(i,4);
                no4= Elements(i,5);
                no5= Elements(i,6);
                no6= Elements(i,7);
                no7= Elements(i,8);
                no8= Elements(i,9);
                edofs=[no1 no5 no2 no6 no3 no7 no4 no8]; % ordem para desenhar    

                fill (coord(edofs,2),coord(edofs,3),pressao(i));hold on
                plot(coord(edofs,2),coord(edofs,3),'r');hold on
                colorbar
            end
            title('Pressão');
            hold off
        case 6 % Sair
            disp('Exiting the program.');
            break;
        case 4 % Erro entre o potencial e o matlab
            figure;
             for i=1:length(Elements)
                no1= Elements(i,2);
                no2= Elements(i,3);
                no3= Elements(i,4);
                no4= Elements(i,5);
                no5= Elements(i,6);
                no6= Elements(i,7);
                no7= Elements(i,8);
                no8= Elements(i,9);
                edofs=[no1 no5 no2 no6 no3 no7 no4 no8]; % ordem para desenhar    

                fill (coord(edofs,2),coord(edofs,3),EPotencial(edofs));hold on
                plot(coord(edofs,2),coord(edofs,3),'r');hold on
                colorbar
             end
            title('Erro do Potencial');
            plot(coord(1),coord(2),'ro');
            hold off
        case 5 % Erro entre a Velocidade
            figure;
             for i=1:length(Elements)
                no1= Elements(i,2);
                no2= Elements(i,3);
                no3= Elements(i,4);
                no4= Elements(i,5);
                no5= Elements(i,6);
                no6= Elements(i,7);
                no7= Elements(i,8);
                no8= Elements(i,9);
                edofs=[no1 no5 no2 no6 no3 no7 no4 no8]; % ordem para desenhar    

                fill (coord(edofs,2),coord(edofs,3),EVelocidade(edofs));hold on
                plot(coord(edofs,2),coord(edofs,3),'r');hold on
                colorbar
             end
            title('Erro da Velocidade');
            plot(coord(1),coord(2),'ro');
            hold off
        otherwise
            disp('Inválido');
    end
end