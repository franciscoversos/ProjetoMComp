%%-------------------------------------------------------------------------
% Flow problem simulation using linear quadrilateral (Q4) elements
%%-------------------------------------------------------------------------
% Import the Nx files that contain the nodes and elements
format long
clear all
fileID0 = fopen('geral_quad4_nodes.txt', 'r');
fileID1 = fopen('geral_quad4_elements.txt', 'r');
Q4NxPotencial  = xlsread("Q4_malhabase_potencial.xlsx");
Q4NxVelocidade = xlsread("Q4_malhabase_velocidade.xlsx");
formatSpec = '%c';'%d' ;
%%-------------------------------------------------------------------------
% Extract the node coordinates, stored in the vector "coord"
nodes = splitlines(fscanf(fileID0, formatSpec));
auxN =[];
for l =11:6:length(nodes)
    auxN = [auxN, nodes(l)];
end
auxN = auxN';
col = split(auxN);
El1=[];
El2 =[];
coord=[];
for l =1:length(col)
    El1 = [El1; str2double(col(l,5))];
    El2 = [El2; str2double(col(l,6))];
    coord = [coord; l, El1(l), El2(l)];
end
%%-------------------------------------------------------------------------
% Build the connectivity matrix: for each element (column 1),
% store its corresponding nodes in the "Elements" array
Big = splitlines(fscanf(fileID1, formatSpec));
auxE =[];
dauxE=[];
for l = 17:length(Big)-4
    auxE = [auxE,Big(l)];
end
auxE = auxE';
colE = split(auxE);
El1 =[];
El2=[];
El3 =[];
El4=[];
Elements=[];
for l =1:length(colE)
    El1 = [El1; str2double(colE(l,11))];
    El2 = [El2; str2double(colE(l,12))];
    El3 = [El3; str2double(colE(l,13))];
    El4 = [El4; str2double(colE(l,14))];
    Elements = [Elements; l, El1(l), El2(l), El3(l),El4(l)];
end

Kr =sparse(length(coord), length(coord));
fr =sparse(length(coord),1);
%%--------------------------------------------------------------------------
% Assemble stiffness matrix (Kr) and force vector (fr)
for l =1:length(Elements)
    no1= Elements(l,2);
    no2= Elements(l,3);
    no3= Elements(l,4);
    no4= Elements(l,5);

    edofs = [no1 no2 no3 no4];
    XN = [coord(no1,2), coord(no1,3);coord(no2,2), coord(no2,3);coord(no3,2), coord(no3,3);coord(no4,2), coord(no4,3)];
    f1 =0;
    [Ke fe B] = Elem_Quad4(XN, f1);
    Kr(edofs,edofs)=Kr(edofs,edofs)+Ke;
    fr(edofs,1)=fr(edofs,1)+fe;
end
%%-------------------------------------------------------------------------
% Boundary conditions:
% Create boundary vectors:
%   - fronteira2: top wall
%   - fronteira3: inner hole
k = boundary(coord(:,2), coord(:,3), 0.6);
fronteira =[coord(k,1), coord(k,2), coord(k,3)];
fronteira(length(fronteira),:)=[];
plot(fronteira(:,2),fronteira(:,3))
coordfront= [];
fronteira1=[];
fronteira2=[];
fronteira3=[];
fronteira4=[];
for l =1:length(fronteira)
    if fronteira(l,2)==0
        fronteira1 = [fronteira1;l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
end
for l =1:length(fronteira)
    if fronteira(l,3)==500
        fronteira2=[fronteira2; l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
    if fronteira(l,2)>=1000 && fronteira(l,2)<=1600 &&fronteira(l,3)>=200 &&fronteira(l,3)<500
        fronteira2=[fronteira2;l, fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
end

for l=1:length(fronteira)
    if fronteira(l,3)==0
        fronteira3=[fronteira3;l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
    if fronteira(l,2) >1750 && fronteira(l,2)<2000 &&fronteira(l,3)<=125&& fronteira(l,3)>0
        fronteira3=[fronteira3;l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
end

for l =1:length(fronteira)
    if fronteira(l,2)==2000
        fronteira4 = [fronteira4;l,fronteira(l,1),fronteira(l,2),fronteira(l,3)];
    end
end
fronteira1=sortrows(fronteira1,4);
fronteira2=sortrows(fronteira2,3);
fronteira3=sortrows(fronteira3,3);
fronteira4=sortrows(fronteira4,3);

% Dirichlet conditions
for l =1:length(fronteira)
    if fronteira(l,2) == 2000
        coordfront = [coordfront;fronteira(l,1), 0];
    end
end
boom = 10^30;
for l = 1:length(coordfront)
    Kr(coordfront(l,1),coordfront(l,1))= boom;
    fr(coordfront(l,1))= boom*coordfront(l,2);
end

% Neumann conditions
% u = 2500;
n=0;
k=0;
nEl=0;
for l =1:length(fronteira)
    if fronteira(l,2)==0 && (fronteira(l,3) == 500 || fronteira(l,3)==0)
        n=n+1;
        fr(fronteira(l,1)) =  2500/2;
    end
    if fronteira(l,2)==0 && (fronteira(l,3)< 500 && fronteira(l,3)>0)
        k=k+1;
        fr(fronteira(l,1)) =  2500;
    end
    nEl=(k+n-1);
end
fr(fronteira(:,1)) =fr(fronteira(:,1))*(500/nEl);
%%-------------------------------------------------------------------------
% Solve potential using backslash
u = Kr\fr;
%%-------------------------------------------------------------------------
% After exporting the Nx data file, compute the error between the two results
EPotencial=[];
for i =1:length(coord)
    EPotencial=[EPotencial; u(i)-Q4NxPotencial(i,8)*10^6];
end

%%-------------------------------------------------------------------------
% Compute pressure and gradients in x and y:
% x-velocity = um, y-velocity = vm, measured at element centroids (xm, ym)
V=[];
for l =1:length(Elements)
    no1= Elements(l,2);
    no2= Elements(l,3);
    no3= Elements(l,4);
    no4= Elements(l,5);
    x1 =coord(no1,2);
    y1=coord(no1,3);
    x2=coord(no2,2);
    y2=coord(no2,3);
    x3=coord(no3,2);
    y3=coord(no3,3);
    x4=coord(no4,2);
    y4=coord(no4,3);
    Ue=[u(no1);u(no2);u(no3);u(no4)];
    XN = [coord(no1,2), coord(no1,3);coord(no2,2), coord(no2,3);coord(no3,2), coord(no3,3);coord(no4,2), coord(no4,3)];
    f1 =0;
    [Ke fe B psi] = Elem_Quad4(XN, f1);
    B =B';
    %
    xm(l)= (x1+x2+x3+x4)/4;
    ym(l)= (y1+y2+y3+y4)/4;
    %
    um(l) = -B(1, :)*Ue;
    vm(l) = -B(2, :)*Ue;
    V= [V;sqrt(um(l)^2 + vm(l)^2)];
    pressao(l)= (0.5).*(2500^2 -(um(l)^2)+ (vm(l)^2));
end

%%-------------------------------------------------------------------------
% Compute velocity at points of interest (here: nodes 1..8)
fno1 = 1;fno2=2;fno3 =3;fno4=4;
fno5 = 5;fno6=6;fno7 =7;fno8=8;

fno=[fno1 fno2 fno3 fno4 fno5 fno6 fno7 fno8];
fno=fno';
W=0;
K=[];
El=[];
U=zeros(length(fno),1);
for l =1:length(coord)
    W=0;
    ufno=0;
    vfno=0;
    for i = 1:length(Elements)
        for k=2:5
            if coord(l,1) == Elements(i,k)
                W=W+1;
                ufno = ufno + um(1,Elements(i,1));
                vfno = vfno + vm(1,Elements(i,1));
            end
        end
    end
    U(l,:)=[sqrt((ufno/W)^2 + (vfno/W)^2)];
    U(fno);
    u(fno);
end
%%----------------------------------------------------------------------
% Velocity error between MATLAB and Nx
EVelocidade=[];
for i=1:length(coord)
    EVelocidade=[EVelocidade;U(i)-Q4NxVelocidade(i,11)*10^9];
end

%%-------------------------------------------------------------------------
% Resultant force on fronteira1
D1=[];
f1x=[];
f1y=[];
d=0;

for l =1:length(fronteira1)-1
    p1 =[fronteira1(l,3); fronteira1(l,4)];
    p2 =[fronteira1(l+1,3);fronteira1(l+1,4)];
    dh=p2-p1;
    d =norm(p2-p1);
    dv = [-dh(2)/d;dh(1)/d];
    D1=[D1,d];
    fx=0;
    for i=1:length(Elements)
        for k =2:5
            for y=2:5
                if fronteira1(l,2)== Elements(i,k)&&fronteira1(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy =pressao(Elements(i,1))*dv(2);
                end
            end
        end
    end
    f1x=[f1x;fx];
    f1y=[f1y;fy];
end
F1x=sum(f1x);
F1y=sum(f1y);
%%----------------------------------------------------------------------
% Resultant force on fronteira2
D2=[];
f2x=[];
f2y=[];
d=0;
dh=0;
for l =1:length(fronteira2)-1
    p1 =[fronteira2(l,3); fronteira2(l,4)];
    p2 =[fronteira2(l+1,3);fronteira2(l+1,4)];
    dh= p2-p1;
    d =norm(p2-p1);
    dv = [-dh(2)/d;dh(1)/d];
    D2=[D2,d];
    fx=0;
    for i=1:length(Elements)
        for k =2:5
            for y=2:5
                if fronteira2(l,2)== Elements(i,k)&&fronteira2(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy = pressao(Elements(i,1))*dv(2);
                end
            end
        end
    end
    f2x=[f2x;fx*d];
    f2y=[f2y;fy*d];
end
F2x=sum(f2x);
F2y=sum(f2y);
%%----------------------------------------------------------------------
% Resultant force on fronteira3 (inner hole)
D3=[];
f3x=[];
f3y=[];
p1=0;
p2=0;
d=0;

for l =1:length(fronteira3)-1
    p1 =[fronteira3(l,3); fronteira3(l,4)];
    p2 =[fronteira3(l+1,3);fronteira3(l+1,4)];
    dh=p2-p1;
    d =norm(p2-p1);
    dv=[dh(2)/d;-dh(1)/d];      % outward normal vector
    D3=[D3,d];
    fx=0;
    for i=1:length(Elements)
        for k =2:5
            for y=2:5
                if fronteira3(l,2)== Elements(i,k)&&fronteira3(l+1,2)==Elements(i,y)
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
% Resultant force on fronteira4
D4=[];
f4x=[];
f4y=[];
d=0;
fx=0;
fy=0;
for l =1:length(fronteira4)-1
    p1 =[fronteira4(l,3); fronteira4(l,4)];
    p2 =[fronteira4(l+1,3);fronteira4(l+1,4)];
    dh= p2-p1;
    d =norm(p2-p1);
    dv=[dh(2)/d;-dh(1)/d];
    D4=[D4,d];
    fx=0;
    for i=1:length(Elements)
        for k =2:5
            for y=2:5
                if fronteira4(l,2)== Elements(i,k)&&fronteira4(l+1,2)==Elements(i,y)
                    fx=pressao(Elements(i,1))*dv(1);
                    fy =pressao(Elements(i,1))*dv(2);
                end
            end
        end
    end
    f4x=[f4x;fx*d];
    f4y=[f4y;fy*d];
end
F4x=sum(f4x);
F4y=sum(f4y);

% Total resultant force
Fx =(F1x +F2x +F3x+F4x)*10^-6;
Fy =(F1y + F2y +F3y + F4y)*10^-6;
%%-------------------------------------------------------------------------
% Menu to choose which plot the user wants to see
while true
    % Menu
    msg  = 'Choose the plot:';
    opts = {'Velocity field', 'Potential', 'Pressure', 'Potential error', 'Velocity error', 'Exit'};
    choice = menu(msg, opts);

    switch choice
        case 1 % Velocity field
            figure;
            hold on;
            title('Velocity field');
            plot(coord(1), coord(2), 'ro');
            quiver(xm, ym, um, vm, 'k');
            plot(fronteira(:,2), fronteira(:,3),'black')
            hold off;

        case 2 % Potential
            figure;
            hold on;
            title('Potential field');
            for l = 1:length(Elements)
                no1 = Elements(l, 2);
                no2 = Elements(l, 3);
                no3 = Elements(l, 4);
                no4 = Elements(l, 5);

                edofs = [no1 no2 no3 no4];
                fill(coord(edofs, 2), coord(edofs, 3), u(edofs));
                hold on;
                colorbar
                plot(coord(edofs, 2), coord(edofs, 3), 'k');
                hold on;
            end
            plot(coord(1), coord(2), 'ro');
            hold off;

        case 3 % Pressure
            figure;
            hold on;
            title('Pressure field');
            for l = 1:length(Elements)
                no1 = Elements(l, 2);
                no2 = Elements(l, 3);
                no3 = Elements(l, 4);
                no4 = Elements(l, 5);

                edofs = [no1 no2 no3 no4];
                fill(coord(edofs, 2), coord(edofs, 3), pressao(l));
                hold on;
                plot(coord(edofs, 2), coord(edofs, 3), 'k');
                hold on;
            end
            colorbar;
            plot(coord(1), coord(2), 'ro');
            hold off;

        case 4 % MATLAB vs NX potential error
            figure;
            hold on;
            title('Potential error');
            for l = 1:length(Elements)
                no1 = Elements(l, 2);
                no2 = Elements(l, 3);
                no3 = Elements(l, 4);
                no4 = Elements(l, 5);

                edofs = [no1 no2 no3 no4];
                fill(coord(edofs, 2), coord(edofs, 3), EPotencial(edofs));
                hold on;
                colorbar
                plot(coord(edofs, 2), coord(edofs, 3), 'k');
                hold on;
            end
            plot(coord(1), coord(2), 'ro');
            hold off;

        case 5 % MATLAB vs NX velocity error
            figure;
            hold on;
            title('Velocity error');
            for l = 1:length(Elements)
                no1 = Elements(l, 2);
                no2 = Elements(l, 3);
                no3 = Elements(l, 4);
                no4 = Elements(l, 5);

                edofs = [no1 no2 no3 no4];
                fill(coord(edofs, 2), coord(edofs, 3), EVelocidade(edofs));
                hold on;
                colorbar
                plot(coord(edofs, 2), coord(edofs, 3), 'k');
                hold on;
            end
            plot(coord(1), coord(2), 'ro');
            hold off;

        case 6 % Exit
            disp('Exiting the program.');
            break;

        otherwise
            disp('Invalid option.');
    end
end
