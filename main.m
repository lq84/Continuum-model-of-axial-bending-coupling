
% This program is modified from the code of the paper
% "Size effects in lattice structures and a comparison to micropolar
% elasticity" by Marcus Yoder, Lonny Thompson, and Joshua Summers.

clear
close all

%% Set up the parameters.

global L e22 nely nelx D macroWidth macroHeight nn IEN coords

%Calculate the Micropolar continuum properties.
L=90; % Total length, mm
e22=-0.015; % Applied strain
N=8; a=L/N; t=0.1*a; theta=30;
D = MP_Square_Mx(a,t,theta)
% D = MP_Square_Mx_Yuan(a,t,theta)

DeformationFactor = 5;

%% Define all of the output fields
P = defaultRun(1);

%If neccesary, set the number of elements in the X and Y direction to
%the default value of 35.
if ~isfield(P,'nelx')
    P.nelx = 5;
end
if ~isfield(P,'nely')
    P.nely = 5;
end
if isempty(P.nelx)
    P.nelx = 5;
end
if isempty(P.nely)
    P.nely = 5;
end

nelx = P.nelx; nely = P.nely;

macroWidth = P.cellWidth*P.numCellsHorz;
macroHeight = P.cellHeight*P.numCellsVert;

%% Make a mesh
if strcmp(P.probType,'halfspace')
    [coords,IEN] = altHalfspaceMesh();
else
    [coords,IEN] = rectMesh();
end
ne = size(IEN,1);  % Number of elements

%% Initialize the global force, k, and B_stored matrixes and boundary conds
nn = size(coords,1);
ndof = nn*3; % Number of degrees of freedom. 3 per node. 
K = sparse(ndof,ndof);
%Create empty matrices for the energy decomposition.
%KNorm = K;
%KShear = K;
%KCurv = K;

%Create empty matrices for the unit cell decomposition
UCK = cell(P.numCellsHorz*P.numCellsVert,1);
for i = 1:P.numCellsHorz*P.numCellsVert
    UCK{i} = K;
end

if P.probType(1)=='F'
    [u,free,essential,~] = FboundaryConditions(coords,P.probType,contProps(P),false);
else
    [u,free,essential,~] = boundaryConditions(coords,P.probType,false);
end
uPresc = u;
%uPrescribed will be used for sorting out which nodes are on the applied 
%displacement face.
%% Generate the local k and f matrixes and solve
%     Add them to the global matrix
% % loop over the elements
etaRow = [-sqrt(3/5),0,sqrt(3/5)]; etaRow = [etaRow,etaRow,etaRow];
xiRow = [-sqrt(3/5),-sqrt(3/5),-sqrt(3/5),0,0,0,sqrt(3/5),sqrt(3/5),sqrt(3/5)];
weight = [25,40,25,40,64,40,25,40,25]/81;
for e = 1:ne
      % loop over local node numbers to get their node global node numbers
      coord = coords(IEN(e,:),:);
      
      % ----------------------
      % Calculate the element stiffness matrix each time.
      % ----------------------
      
      ke = sparse(8*3,8*3);
      
      % Loop over the guass points
      for gu = 1:9
          eta = etaRow(gu);
          xi = xiRow(gu);
          wght = weight(gu);
          [B, J_det] = BMatrixQ(xi,eta,coord);
          ke = ke + transpose(B)*D*B*J_det*wght;
      end
           
     % Insert the element stiffness matrix into the global.
     nodes1 = IEN(e,:);
     % I cannot use the union, or else the order get messed up. The order
     % is important. Same in the actual topology code when you are
     % calculating the objective.
     dofNumbers(1:3:24) = nodes1*3-2;
     dofNumbers(2:3:24) = nodes1*3-1;
     dofNumbers(3:3:24) = nodes1*3-0;
     % multiply t*ke and add to global matrix. t = thickness or t_z
     K(dofNumbers,dofNumbers) = K(dofNumbers,dofNumbers) + ke;%#ok
     
     %Repeat for the unit cell matrices.
     whichUC = floor(mean(coords(nodes1,:))./[P.cellWidth,P.cellHeight])+1;
     whichUC = (whichUC(2)-1)*P.numCellsHorz+whichUC(1);
     UCK{whichUC}(dofNumbers,dofNumbers) = UCK{whichUC}(dofNumbers,dofNumbers) + ke;
end

K_ff = K(free,free);
F_f = zeros(length(free),1);
u(free) = K_ff \ (F_f-K(free,essential)*u(essential));

P.strainEnergyCont = u'*K*u/2;



if strcmp(P.probType,'halfspace')
    DeformationFactor = -2;
end

%% Do the plots.
if or(strcmp(P.plots,'all'),not(isempty(strfind(P.plots,'displacement'))))
    figure
    hold on
    sortingVec = [1,5,2,6,3,7,4,8];
    %This turns the element numbered coordinates into ccw numbered coordinates.
    %Plot the deformed configuration
    for e = 1:ne
        % loop over local node numbers to get their node global node numbers
        % get the global X,Y position of each node and put in array
        coord = coords(IEN(e,sortingVec),:);
        
        %% Plot the element outline and the displacments
        coordD = zeros(9,2);
        nodes  = IEN(e,sortingVec);
        for temp = 1:8
%             text(coord(temp,1),coord(temp,2),int2str(nodes(temp))); % Node number
            coordD(temp,2) = coord(temp,2) + DeformationFactor*u(3*nodes(temp)-1); % Y value
            coordD(temp,1) = coord(temp,1) + DeformationFactor*u(3*nodes(temp)-2); % X value 
            %coordD(temp,1) =  coordD(temp,1) - coordD(temp,2)*.01*multiplierScale; %This line deshears the shear.
        end
        
        coordI = coord;
        coordD(9,:) = coordD(1,:);
        coordI(9,:) = coordI(1,:);
        plot(coordI(:,1),coordI(:,2),'-g');
        plot(coordD(:,1),coordD(:,2), '-b');
    end
    %Make axis settings for the element outlines.
    if macroWidth>macroHeight
        xlim([-0.0075*DeformationFactor,1+0.0125*DeformationFactor]*macroWidth)
    else
        ylim([-0.0075*DeformationFactor,1+0.0125*DeformationFactor]*macroHeight)
    end
    axis equal
    tti= strcat('Element Deformation - Displacement of the elements shown in Blue');
    title(tti);
    hold off
end

figure
hold on
sortingVec = [1,5,2,6,3,7,4,8];
maxU1=0;
maxU2=0;
%This turns the element numbered coordinates into ccw numbered coordinates.
%Plot the deformed configuration
for e = 1:ne
    % loop over local node numbers to get their node global node numbers
    % get the global X,Y position of each node and put in array
    coord = coords(IEN(e,sortingVec),:);
    
    %% Plot the element outline and the displacments
    coordD = zeros(9,2);
    dispy=zeros(9,1); dispx=zeros(9,1); disp=zeros(9,1); 
    nodes  = IEN(e,sortingVec);
    for temp = 1:8
        coordD(temp,2) = coord(temp,2) + DeformationFactor*u(3*nodes(temp)-1); % Y value
        coordD(temp,1) = coord(temp,1) + DeformationFactor*u(3*nodes(temp)-2); % X value
        dispy(temp)=u(3*nodes(temp)-1);
        dispx(temp)=u(3*nodes(temp)-2);
        
        disp(temp)=sqrt(dispy(temp)^2+dispx(temp)^2);
        %coordD(temp,1) =  coordD(temp,1) - coordD(temp,2)*.01*multiplierScale; %This line deshears the shear.
    end
    
    coordD(9,:) = coordD(1,:);
    dispy(9)=dispy(1); dispx(9)=dispx(1); disp(9)=disp(1); 
    maxU1=max(maxU1,max(abs(dispx)));
    maxU2=max(maxU2,max(abs(dispy)));
    patch(coordD(:,1),coordD(:,2),dispx,'EdgeColor','none');
end
%Make axis settings for the element outlines.
if macroWidth>macroHeight
    xlim([-0.0075*DeformationFactor,1+0.0125*DeformationFactor]*macroWidth)
else
    ylim([-0.0075*DeformationFactor,1+0.0125*DeformationFactor]*macroHeight)
end
axis equal
tti= strcat('Displacement Contour, u1');
title(tti);
colormap jet
colorbar
hold off
axis off
% k=8*maxU1/L^2   %unit: mm^-2
ka_over_e=8*maxU1/N/maxU2

nodes_vertical=9:26:217;
nodes_ycoord=coords(nodes_vertical,2);
nodes_u1=u(nodes_vertical*3-2);
P= polyfit(nodes_ycoord, nodes_u1, 2);
k=2*P(1)
xi=linspace(min(nodes_ycoord),max(nodes_ycoord),20);  
yi= polyval(P, xi);
figure
plot(nodes_ycoord,nodes_u1,'b*',xi,yi,'k');
view(90,-90)
legend('Continuum model','Curve fitting','Location','east')
RFnodes=K*u;
RF2=sum(abs(RFnodes(2:3:end)))/2





function [coords,IEN] = rectMesh()
global nelx nely macroWidth macroHeight nn
% This corresponds to an automatic rectangular mesh.
    %ne = nelx*nely; % number of elements
    nn = (2*nelx+1)*(nely+1)+nely*(nelx+1); % number of nodes

    el = 1;%In this loop count stores the element number.
    nnRow = 2*nelx+1;
    nnHalfRow = nelx+1;
    nnFullRow = nnRow+nnHalfRow;
    %Node number mapping stsarts with bottom left corner and goes through
    %the row, then the half row, then the next row. The first node in the
    %i-th full row is (i-1)*(nnRow+nnHalfRow)+1
    %nnCol = nely+1;
    IEN = zeros(nelx*nely,8); % Index of element nodes (IEN)
    % Each row, so nely # of row
    for i = 1:nely
        zerothNode = (i-1)*(nnRow+nnHalfRow);
        fullRowEnd = zerothNode+nnRow;
        % Each column, so nelx # of row
        for j= 1:nelx
            IEN(el,[1,5,2])=zerothNode+(j-1)*2+(1:3);
            IEN(el,[8,6])= fullRowEnd + j + [0,1];
            IEN(el,[4,7,3])=zerothNode+(j-1)*2+(1:3)+nnFullRow;
            el = el+1;
        end
    end

    % Find and store the global positions of each node
    % Each element rectangle. The aspect ratio is determined by the 
    % number off elements in the X and Y directions and the height and 
    % width of the beam.
    %
    % Store both the X and Y positions
    coords = zeros(nn,2);
    nextNode = 1; 
    for row = 1:nely
        %Enter the coordinates for the row.
        coords(nextNode+(1:nnRow)-1,1) = linspace(0,macroWidth,nnRow);
        coords(nextNode+(1:nnRow)-1,2) = (row-1)*macroHeight/nely;
        nextNode = nextNode+nnRow;
        %Enter the coordinates for the half row.
        coords(nextNode+(1:nnHalfRow)-1,1) = linspace(0,macroWidth,nnHalfRow);
        coords(nextNode+(1:nnHalfRow)-1,2) = (row-.5)*macroHeight/nely;
        nextNode = nextNode+nnHalfRow;
        %plot(coords(:,1),coords(:,2),'o'); xlim([-1,macroWidth+1]); ylim([-1,macroHeight+1]); 
    end
    %Enter the coordinates for the top row.
    coords(nextNode+(1:nnRow)-1,1) = linspace(0,macroWidth,nnRow);
    coords(nextNode+(1:nnRow)-1,2) = macroHeight;
end