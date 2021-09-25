function P = continuumSolverQ(varargin)
%This function takes in a parameters structure, P, and solves the
%appropriate micropolar FEA. The output structure has all of the relevant
%postprocessed results from the latticeSolver filled in.
%If there are no input arguments, it will run whatever problem is defined in
%the defaultRun function.

%tic;

%% Set up the parameters.
%This way the solver can be called multiple times with slightly different
%parameters and it.

global nely nelx D macroWidth macroHeight nn IEN coords
if isempty(varargin)
    %Define all of the output fields
    P = defaultRun(1);
else
    P = varargin{1};
end

%Calculate the continuum properties.
D = contProps(P);

%If neccesary, set the number of elements in the X and Y direction to
%the default value of 35.
if ~isfield(P,'nelx')
    P.nelx = 35;
end
if ~isfield(P,'nely')
    P.nely = 35;
end
if isempty(P.nelx)
    P.nelx = 35;
end
if isempty(P.nely)
    P.nely = 35;
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
ne = size(IEN,1);

%% Initialize the global force, k, and B_stored matrixes and boundary conds
nn = size(coords,1);
ndof = nn*3; % Number of degrees of freedom. 2 per node. 
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

multiplierScale = 10;
if strcmp(P.probType,'halfspace')
    multiplierScale = -2;
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
            coordD(temp,2) = coord(temp,2) + multiplierScale*u(3*nodes(temp)-1); % Y value
            coordD(temp,1) = coord(temp,1) + multiplierScale*u(3*nodes(temp)-2); % X value 
            %coordD(temp,1) =  coordD(temp,1) - coordD(temp,2)*.01*multiplierScale; %This line deshears the shear.
        end
        
        %coord2 = coord;
        coordD(9,:) = coordD(1,:);
        %coord2(9,:) = coord2(1,:);
        %plot(coord2(:,1),coord2(:,2),'-g');
        plot(coordD(:,1),coordD(:,2), '-b');
    end
    %Make axis settings for the element outlines.
    if macroWidth>macroHeight
        xlim([-0.0075*multiplierScale,1+0.0125*multiplierScale]*macroWidth)
    else
        ylim([-0.0075*multiplierScale,1+0.0125*multiplierScale]*macroHeight)
    end
    axis equal
    tti= strcat('Element Deformation - Displacement of the elements shown in Blue');
    title(tti);
    hold off
end

end

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