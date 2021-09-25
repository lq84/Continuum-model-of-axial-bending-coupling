function P = classicalSolverQ(varargin)
%This function takes in a parameters structure, P, and solves the
%appropriate lattice FEA. The output structure has all of the relevant
%postprocessed results from the latticeSolver filled in.
%If there are no input arguments, it will run whatever problem is defined in
%the defaultRun function.
tic;

%% Set up the parameters.
%This way the solver can be called multiple times with slightly different
%parameters and it.


global D macroWidth macroHeight nn IEN coords
if isempty(varargin)
    %Define all of the output fields
    P = defaultRun(1);
else
    P = varargin{1};
end

if length(varargin)==6
    D = varargin{6};
else
    D = contProps(P);
end
%% Skip the boundary conditions curve 1 and curve 2
%this block of code needs to produce results that won't cause other blocks
%of code to crash. 
if strcmp('curve',P.probType(1:5))
    P.strainEnergyClass = 0;
    return
end
%Strip the micropolar behavior.
mu = (D(3,3)+D(3,4))/2;
D(3:6,1:6) = 0; D(1:6,3:6) = 0;
D(3,3) = mu;
D = D(1:3,1:3);

macroWidth = P.cellWidth*P.numCellsHorz; 
macroHeight = P.cellHeight*P.numCellsVert;

%% Make a mesh
%First check that the number of elements has been defined.
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

if strcmp(P.probType,'halfspace')
    [coords,IEN] = altHalfspaceMesh();
else
    %Expand the number of elements to make sure the energy decompostion is 
    %the same.
    %P.nelx = ceil(P.nelx/P.numCellsHorz)*P.numCellsHorz;
    %P.nely = ceil(P.nely/P.numCellsVert)*P.numCellsVert;
    [coords,IEN] = rectMesh(P.nelx,P.nely);
end
ne = size(IEN,1);
%% Initialize the global force, k, and B_stored matrixes and boundary conds
nn = size(coords,1);
ndof = nn*2; % Number of degrees of freedom. 2 per node. 
K = sparse(ndof,ndof);

%Create empty matrices for the unit cell decomposition
UCK = cell(P.numCellsHorz*P.numCellsVert,1);
for i = 1:P.numCellsHorz*P.numCellsVert
    UCK{i} = K;
end

[u,free,essential,~] = boundaryConditions(coords,P.probType,false);
%Strip the micropolar DoFs from these vectors.
nonMPDofs = true(nn*3,1);
nonMPDofs(3:3:end) = false;
uPresc = u(nonMPDofs,1);
u = uPresc;

%Strip micropolar dofs and Change the dof numbers in free
free = free(not(mod(free,3)==0)); %strip micropolar dofs
xDofs = mod(free,3)==1; yDofs = mod(free,3)==2; %mark the nodes that need changing.
free(xDofs) = 2*(free(xDofs)+2)/3-1;
free(yDofs) = 2*(free(yDofs)+1)/3;
%Strip micropolar dofs and Change the dof numbers in essential
essential = essential(not(mod(essential,3)==0));
xDofs = mod(essential,3)==1; yDofs = mod(essential,3)==2;
essential(xDofs) = 2*(essential(xDofs)+2)/3-1;
essential(yDofs) = 2*(essential(yDofs)+1)/3;

%uPrescribed will be used for sorting out which nodes are on the applied 
%displacement face.
%% Generate the local k and f matrixes and solve
% Add them to the global matrix
% loop over the elements
etaRow = [-sqrt(3/5),0,sqrt(3/5)]; etaRow = [etaRow,etaRow,etaRow];
xiRow = [-sqrt(3/5),-sqrt(3/5),-sqrt(3/5),0,0,0,sqrt(3/5),sqrt(3/5),sqrt(3/5)];
weight = [25,40,25,40,64,40,25,40,25]/81;
for e = 1:ne
      % loop over local node numbers to get their node global node numbers
      coord = coords(IEN(e,:),:);
      
      % Calculate the element stiffness matrix each time.      
      ke = sparse(8*2,8*2);
           
      % Loop over the guass points
      for gu = 1:9
          eta = etaRow(gu);
          xi = xiRow(gu);
          wght = weight(gu);
          [B, J_det] = BClass(xi,eta,coord);
          ke = ke + transpose(B)*D*B*J_det*wght;
      end
     
     % Insert the element stiffness matrix into the global.
     nodes1 = IEN(e,:);
     % I cannot use the union, or else the order get messed up. The order
     % is important. Same in the actual topology code when you are
     % calculating the objective.
     dofNumbers(1:2:16) = nodes1*2-1;
     dofNumbers(2:2:16) = nodes1*2-0;
     
     % multiply t*ke and add to global matrix. t = thickness or t_z
     K(dofNumbers,dofNumbers) = K(dofNumbers,dofNumbers) + ke;%#ok
     
     %Repeat for the unit cell matrices.
     whichUC = floor(mean(coords(nodes1,:))./[P.cellWidth,P.cellHeight])+1;
     whichUC = (whichUC(2)-1)*P.numCellsHorz+whichUC(1);
     UCK{whichUC}(dofNumbers,dofNumbers) = UCK{whichUC}(dofNumbers,dofNumbers) + ke;
end

% %The micropolar DoFs are no longer free or essential.
% free = free(not(mod(free,3)==0));
% essential = essential(not(mod(essential,3)==0));

K_ff = K(free,free);
F_f = zeros(length(free),1);
u(free) = K_ff \ (F_f-K(free,essential)*u(essential));

%% Calculate the strain energy per unit cell.
% %Find out which nodes are in which unit cells.
% if strcmp('halfspace',P.probType)
%     P.UCStrainEnergyClass = ones(1,P.numCellsHorz*P.numCellsVert);
% else
%     P.UCStrainEnergyClass = zeros(1,P.numCellsHorz*P.numCellsVert);
%     for i = 1:P.numCellsHorz*P.numCellsVert
%         P.UCStrainEnergyClass(i) = u'*UCK{i}*u/2;
%     end
% end
%% Use the energy decomposition stiffness matrices.
P.strainEnergyClass = 1/2*u'*K*u;

multiplierScale = 10;
if strcmp(P.probType,'halfspace')
    multiplierScale = 0;
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
            coordD(temp,2) = coord(temp,2) + multiplierScale*u(2*nodes(temp)-0); % Y value
            coordD(temp,1) = coord(temp,1) + multiplierScale*u(2*nodes(temp)-1); % X value 
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

function essential = halfspaceMesh() %#ok
%This patch of code will create a mesh appropriate for a micropolar half space.
global nelx L IEN coords nn
%% Set up mesh global geometric parameters.
%The mesh is divided into two regions, a skew grid close to the origin, and
%a polar region farther away from it.
elPerEdgeSG = nelx;
edgeLengthSG = L;
%outsideRad = edgeLengthSG*10;
outsideRad = 200; %This matches the abaqus simulation.
polarBias = 2;
%This describes how big each cell in the polar region is compared to the one before. 

skewness = 1.065;
%This describes how crooked/skew the grid near the origin is. The skewness
%smooths the transition of the mesh the squarish internal mesh and the
%circlish polar region. A skewness of 1 makes all of the 
%corners of the skew grid equally far from the origin. A skewness of
%sqrt(2) makes the grid perfectly square.
%% Set up a skew grid patch of dense mesh close to the origin.
%Set up seed nodes.
%botEdgeSeeds = [linspace(0,edgeLengthSG,elPerEdgeSG+1)',zeros(elPerEdgeSG+1,1)];
leftEdgeSeeds = [zeros(elPerEdgeSG+1,1),linspace(0,edgeLengthSG,elPerEdgeSG+1)'];
farCorner = [sqrt(.5)*skewness*edgeLengthSG , sqrt(.5)*skewness*edgeLengthSG];
rightEdgeSeeds = [linspace(edgeLengthSG,farCorner(1),elPerEdgeSG+1)',...
                  linspace(0,farCorner(2),elPerEdgeSG+1)'];
topEdgeSeeds = [linspace(0,farCorner(1),elPerEdgeSG+1)',...
                linspace(edgeLengthSG,farCorner(2),elPerEdgeSG+1)'];

%Initialize a global position matrix. Before we do that we need to figure
%out how many nodes ought to be in the polar region. In order to do that,
%we need to figure out how big the first ring should be. The radius of the
%first ring should make the element right outside the corner of the skew
%grid have two equal length sides.
dist = (topEdgeSeeds(end-1,:)-farCorner);
dist = sqrt(dist*dist');
insideRad = edgeLengthSG*skewness+dist;

%This next block figures out how many rings are required to make the polar
%space with the proper polarBias.
radNumber = 1; %RadNumber is the number of node radii that are 
radii = insideRad;
while radii < outsideRad
    radii = [radii; radii(end) + radii(end)*pi/2/(elPerEdgeSG*2+1)*polarBias];%#ok grow
    radNumber = radNumber+1;
end
%calculate the number of nodes and initialize the coords matrix.
nn = (elPerEdgeSG+1)^2+(elPerEdgeSG*2+1)*radNumber; %Number of nodes in the global mesh.
coords = zeros(nn,2);
%Initialize the essential dofs vector.
%3 essential dofs per each node on the outside (elPerEdgeSG*2+1)*3
%2 essential dofs per node on the midline, phi and u1 (elPerEdgeSG+1+radNumber)*2
essential = zeros((elPerEdgeSG*2+1)*3 + (elPerEdgeSG+1+radNumber)*2,1);
essentialCount = 1;
%Essential count tells where in the essential vector to put the next essential nodes.
%Interpolate between the nodes.
for i = 1:elPerEdgeSG+1
    temp = [linspace(leftEdgeSeeds(i,1),rightEdgeSeeds(i,1),elPerEdgeSG+1)',...
            linspace(leftEdgeSeeds(i,2),rightEdgeSeeds(i,2),elPerEdgeSG+1)'];
    nodes = (i-1)*(elPerEdgeSG+1)+1:i*(elPerEdgeSG+1);
    coords(nodes,1:2) = temp;
    %Add the nodes on the midline to the essential vector.
    essential([essentialCount,essentialCount+1]) = [nodes(1)*3-2,nodes(1)*3];
    essentialCount = essentialCount+2;
end
%Create the IEN for this region.
IEN = zeros(elPerEdgeSG^2+(radNumber)*elPerEdgeSG*2,4);
count = 1;
for i = 1:elPerEdgeSG
    rowMultiplier = i-1;
    % Each column, so nelx # of row
    for j= 1:elPerEdgeSG
        IEN(count,:)=[rowMultiplier    *(elPerEdgeSG+1)+j,   ...
                      rowMultiplier    *(elPerEdgeSG+1)+j+1, ...
                     (rowMultiplier +1)*(elPerEdgeSG+1)+j+1, ...
                     (rowMultiplier +1)*(elPerEdgeSG+1)+j];
        count = count+1;
    end
end
%% Set up the polar circle patch.
%Start with the nodalPositions.
currentNode = (elPerEdgeSG+1)^2+1;
for theta = linspace(0,pi/2,elPerEdgeSG*2+1)
    lastNode = currentNode+radNumber-1;
    coords(currentNode:lastNode,:) = ...
        [radii*cos(theta),radii*sin(theta)];
    %Add the last node's dofs to the essential vector.
    essential(essentialCount:essentialCount+2) = lastNode*3-2:lastNode*3;
    currentNode = currentNode + radNumber;
    essentialCount = essentialCount+3;
end
%Add the remainder of the midline dofs to the essential vector.
essential(essentialCount:essentialCount+radNumber) = ...
    (nn-radNumber)*3-2:3:nn*3;
essentialCount = essentialCount+radNumber;
essential(essentialCount:essentialCount+radNumber) = (nn-radNumber)*3:3:nn*3;
%Create an IEN for this region.
%Begin by creating a vector with the node numbers that will be linked to 
%from the skew grid region.
danglingNodes = [elPerEdgeSG+1:elPerEdgeSG+1:(elPerEdgeSG+1)^2,...
    ((elPerEdgeSG+1)^2-1):-1:(elPerEdgeSG+1)^2-elPerEdgeSG];
for i = 1:elPerEdgeSG*2
    sgNodes = (elPerEdgeSG+1)^2;
    rowMultiplier = i-1;
    % Each column, so nelx # of row
    for j= 1:radNumber
        if j==1
            IEN(count,:)=[danglingNodes(i),...
                 sgNodes+(rowMultiplier  )*(radNumber)+j, ...
                 sgNodes+(rowMultiplier+1)*(radNumber)+j, ...
                          danglingNodes(i+1)];
        else
            IEN(count,:)=[rowMultiplier    *(radNumber)+j-1,   ...
                         (rowMultiplier   )*(radNumber)+j, ...
                         (rowMultiplier +1)*(radNumber)+j, ...
                         (rowMultiplier +1)*(radNumber)+j-1]+sgNodes;
        end
        count = count+1;
    end
end
%Round any tiny elements of coords to zero.
coords(abs(coords)<10^-6) = 0;
%remove duplicates from essential.
essential = unique(essential);
end
function [coords,IEN] = rectMesh(nelx,nely)
global macroWidth macroHeight nn
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
function [B, J_det] = BClass(xi,eta,coord)
%This sets up the matrix that helps transform the displacement into
%the strain and then the stress.
%The inputs are:
%xi and eta, scalar values.
%coord a 8x2 matrix containing the (x,y) locations of the element's nodes.
%The outputs are the B matrix, a 6x24 matrix, and J_det the determinant of
%the Jacobian.
B_hat = [-1/4*(1-eta)*(-xi-eta-1)-1/4*(1-xi)*(1-eta),...
         -1/4*(1-xi )*(-xi-eta-1)-1/4*(1-xi)*(1-eta);
          1/4*(1-eta)*(xi-eta-1 )+1/4*(1+xi)*(1-eta),...
         -1/4*(1+xi)*(xi-eta-1)-1/4*(1+xi)*(1-eta);
          1/4*(1+eta)*(xi+eta-1)+1/4*(1+xi)*(1+eta),...
          1/4*(1+xi)*(xi+eta-1)+1/4*(1+xi)*(1+eta);
         -1/4*(1+eta)*(-xi+eta-1)-1/4*(1-xi)*(1+eta),...
          1/4*(1-xi )*(-xi+eta-1)+1/4*(1-xi)*(1+eta);
          -xi*(1-eta),-1/2+1/2*xi^2;
          1/2-1/2*eta^2,-eta*(1+xi);
          -xi*(1+eta),1/2-1/2*xi^2;
          -1/2+1/2*eta^2,-eta*(1-xi)]';

% Calculate the Jacobian
J=B_hat*coord;

% Calculate the determinate
J_det = det(J);
% J_inverse = inv(J);

% Form the B matrix
% B_2by4 = J_inverse*B_hat;
%dbstop if warning
B_2by8 = J\B_hat;

% Form B, which is an 6 by 24
% The form of this matrix was derived from
%Quadrilateral isoparametric finite elements for plane elastic
%Cosserat bodies
%My notes on this are written up in a Maple document of the same
%name.
B = sparse(3,8*2);
%corresponds to epsilon_xx = ux,x
B(1,1:2:16) = B_2by8(1,:);

%corresponds to epsilon_yy = uy,y
B(2,2:2:16) = B_2by8(2,:);

%Corresponds to epsilon_xy = uy,x
B(3,2:2:16) = B_2by8(1,:);
B(3,1:2:16) = B_2by8(2,:);
end