function [u,free,essential,F] = boundaryConditions(varargin)
%This function applies takes in a set of nodal coordinates and returns the
%appropriate variables to create boundary conditions.

%% Find the free and essential boundary conditions.
%Notes about notation. There are three dofs per node that can be fixed. 3
%variables are set up for xFixed, yFixed and thetaFixed. These are then
%combined into an essential vector. There is also a vector for xMoves,
%yMoves, thetaMoves as well as three variables for the magnitude of those
%movements.
globalPos = varargin{1};
probType = varargin{2};
isLattice = varargin{3};
global e22

nn = size(globalPos,1);
macroWidth = max(globalPos(:,1));
macroHeight = max(globalPos(:,2));
xFixed = zeros(nn,1);
yFixed = zeros(nn,1);
thetaFixed = zeros(nn,1);
xMoves = zeros(nn,1);     xMag = 0;
yMoves = zeros(nn,1);     yMag = 0;
thetaMoves = zeros(nn,1); thetaMag = 0.01; %0.01 radians
F = zeros(nn*3,1);

uRot = zeros(nn,3);
tol = 10^-5; %Nodes that are this close to the boundary or closer 
        %will be considered to be on the boundary.
switch probType
    case 'stretchX'
        %stretchX holds the nodes on the left end fixed, prevents the
        %nodes on the right end from rotating or moving vertically and
        %applies a 1% displacement to the nodes on the right end. Get
        %vectors of nodes on the right and left.
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
                         globalPos(:,1) >= macroWidth-tol);
        thetaFixed = or(leftNodes,rightNodes);             
        xFixed = leftNodes;
        yFixed = or(leftNodes,rightNodes);
        xMag = macroWidth/100;
        xMoves = rightNodes;
    case 'stretchXfree'
        %This boundary condition will not fix the contraction at the
        %boundaries.
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
                         globalPos(:,1) >= macroWidth-tol);
        %thetaFixed = or(leftNodes,rightNodes);             
        xFixed = leftNodes;
        
        yFixed = and(globalPos(:,1) <= tol,...
                     globalPos(:,2) <= tol);
        xMag = macroWidth/100;
        xMoves = rightNodes;
        
    case 'stretchY'
        %stretchY holds the nodes on the bot end fixed, prevents the
        %nodes on the top end from moving vertically and
        %applies a 1% displacement to the nodes on the right end. Get
        %vectors of nodes on the right and left.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
%         thetaFixed = or(botNodes,topNodes);
        xFixed = or(botNodes,topNodes);
        yFixed = or(botNodes,topNodes);
%         thetaFixed = or(botNodes,topNodes);
        yMoves = topNodes;
        yMag = macroWidth*e22; 
    case 'stretchYfree'
        %stretchYfree is like stretchY but it does not constrain
        %contraction.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
        %thetaFixed = or(botNodes,topNodes);
        xFixed = and(globalPos(:,1) <= tol,...
                     globalPos(:,2) <= tol);
        yFixed = or(botNodes,topNodes);
        %thetaFixed = or(botNodes,topNodes);
        yMoves = topNodes;
        yMag = macroWidth*e22;    
    case 'transverseY'
        %transverseY holds the nodes on the bot end fixed, prevents the
        %nodes on the top end from rotating or moving vertically and
        %applies a 1% x direction displacement to the nodes on the top. 
        
        %Get vectors of nodes on the right and left.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
        thetaFixed = or(botNodes,topNodes);
        xFixed = botNodes;
        yFixed = or(botNodes,topNodes);
        xMoves = topNodes;
        xMag = macroHeight/100;
        
    case 'transverseX'
        %transverse12 holds the nodes on the left end fixed, prevents the
        %nodes on the right end from rotating or moving vertically and
        %applies a 1% y direction displacement to the nodes on the right. 
        
        %Get vectors of nodes on the right and left.
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
                         globalPos(:,1) >= macroWidth-tol);
        thetaFixed = or(leftNodes,rightNodes);
        yFixed = leftNodes;
        xFixed = or(rightNodes,leftNodes);
        thetaFixed = or(leftNodes,rightNodes);
        yMoves = rightNodes;
        yMag = macroWidth/100;
        
    case 'transverseYhinge'
        %transverseY holds the nodes on the bot end fixed, prevents the
        %nodes on the top end from rotating or moving vertically and
        %applies a 1% x direction displacement to the nodes on the top. 
        
        %Get vectors of nodes on the right and left.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
        %thetaFixed = or(botNodes,topNodes);
        xFixed = botNodes;
        yFixed = or(botNodes,topNodes);
        thetaFixed = zeros(nn,1);
        xMoves = topNodes;
        xMag = macroHeight/100;
        
    case 'curveX'
        %Curve1 applies a k13 curvature to the mesh. The left end is fixed.
        %The right end has a rotation applied to it equal to
        %macroWidth*0.01. The right end has the translation dofs fixed.
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
                         globalPos(:,1) >= macroWidth-tol);
        thetaFixed = leftNodes;
        thetaMoves = rightNodes;
        yFixed = or(rightNodes,leftNodes);
        xFixed = or(rightNodes,leftNodes);
        
    case 'curveY'
        %Curve1 applies a k23 curvature to the mesh. The bot end is fixed.
        %The top end has a rotation applied to it equal to
        %macroHeight*0.01. The top end has the translation dofs fixed.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
        thetaFixed = or(botNodes,topNodes);
        thetaMoves = topNodes;
        yFixed = or(botNodes,topNodes);
        xFixed = or(botNodes,topNodes);
        
    case 'halfspace'
        %halfspace holds the nodes on the left side symmetric ux=0 and
        %phi=0, and holds the nodes on the top and right sides fixed. It
        %applies a displacement of .1 to the bottom left corner.
        
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
                         globalPos(:,1) >= macroWidth-tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);

        thetaFixed = or(leftNodes,topNodes);
        xFixed = or(or(leftNodes,rightNodes),topNodes);
        yFixed = or(topNodes,rightNodes);
        yMag = macroHeight/100;
        yMoves(1) = 1; %These vectors are node numbers not dof numbers.
        
    case 'bendYfixed'
        %This boundary condition is like bendY, but the whole bottom side is fixed in the X direction.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
        xFixed = botNodes;
        yFixed = botNodes;
        thetaFixed = botNodes;
        xMean = (max(globalPos(topNodes,1))-min(globalPos(topNodes,1)))/2;
        yMag = 0.01;
        for nodeI=find(topNodes)
            thetaMoves(nodeI) = -1;
            yMoves(nodeI)  = xMean-globalPos(nodeI,1);
        end
    case 'bendXfixed'
        %This boundary condition is like bendX, but the whole left side is fixed in the Y direction.
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
            globalPos(:,1) >= macroWidth-tol);
        yFixed = leftNodes;
        xFixed = leftNodes;
        thetaFixed = leftNodes;
        yMean = (max(globalPos(rightNodes,2))+min(globalPos(rightNodes,2)))/2;
        xMag = -0.01;
        for nodeI=find(rightNodes)
            thetaMoves(nodeI) = -1;
            xMoves(nodeI)  = yMean-globalPos(nodeI,2);
        end
    case 'bendY'
        %This boundary condition is like bendYfree, but only the bottom
        %left corner is fixed in the x direction.
        botNodes = and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol);
        topNodes = and(globalPos(:,2) <= macroHeight+tol,...
                         globalPos(:,2) >= macroHeight-tol);
        xFixed = and(botNodes,and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol));
        yFixed = botNodes;
        thetaFixed = botNodes;
        xMean = (max(globalPos(topNodes,1))+min(globalPos(topNodes,1)))/2;
        yMag = 0.01;
        for nodeI=find(topNodes)
            thetaMoves(nodeI) = -1;
            yMoves(nodeI)  = xMean-globalPos(nodeI,1);
        end
    case 'bendX'
        %This boundary condition is like bendYfree, but only the bottom
        %left corner is fixed in the x direction.
        leftNodes = and((globalPos(:,1) <= tol) ,globalPos(:,1) >= -tol);
        rightNodes = and(globalPos(:,1) <= macroWidth+tol,...
                         globalPos(:,1) >= macroWidth-tol);
        yFixed = and(leftNodes,and((globalPos(:,2) <= tol) ,globalPos(:,2) >= -tol));
        xFixed = leftNodes;
        thetaFixed = leftNodes;
        yMean = (max(globalPos(rightNodes,2))+min(globalPos(rightNodes,2)))/2;
        xMag = -0.01;
        for nodeI=find(rightNodes)
            thetaMoves(nodeI) = -1;
            xMoves(nodeI)  = yMean-globalPos(nodeI,2);
        end
    otherwise
        error('Unrecognized problem type')
end

%% Process the boundary conditions and format them to be returned.
%Set up the essential dofs in a vector.
essential = zeros(nn*3,1);
essential(1:3:end,1) = or(xFixed,xMoves);
essential(2:3:end,1) = or(yFixed,yMoves);
essential(3:3:end,1) = or(thetaFixed,thetaMoves);

free = find(1-essential);
essential = find(essential);
%Set up the uEssential vector. In order for a node to have a perscribed
%non-zero movement, it must be in a nMoves vector and be part of essential.
u = zeros(nn*3,1);
u(1:3:end,1) = xMoves*xMag+uRot(:,1);
u(2:3:end,1) = yMoves*yMag+uRot(:,2);
u(3:3:end,1) = thetaMoves*thetaMag;
%checkBC()

function checkBC()
    %This subfunction will plot out the boundary conditions.
    %Begin by drawing out the square that is the boundary of the region.
    figure
    hold on
    xlim([-macroWidth*0.2,macroWidth*1.2])
    ylim([-macroHeight*0.2,macroHeight*1.2])
    plot([0,macroWidth,macroWidth,0,0],[0,0,macroHeight,macroHeight,0],'k')
    symScale = -0.1*max(max(globalPos)); %symScale regulates how big the symbols are.
    %loop over the nodes and plot them with their symbol.
    %The symbols for xFixed and xMoves are a red and blue T oriented along
    %the appropriate axis.
    %The symbol for thetaFixed is a green quarter circle. 
    %The symbol for xMoves and yMoves are red and blue arrows oriented along
    %the appropriate axis and facing the opposite direction from the fixed
    %symbols.
    %The symbol for thetaMoves is a green double quarter circle facing the
    %opposite direction from the theta fixed circle.
    for node = 1:nn
        if xFixed(node)
            xVec = globalPos(node,1) + [0,1,  1,   1]*symScale;
            yVec = globalPos(node,2) + [0,0,0.4,-0.4]*symScale;
            plot(xVec,yVec,'r');
        end
        if yFixed(node)
            xVec = globalPos(node,1) + [0,0,0.4,-0.4]*symScale;
            yVec = globalPos(node,2) + [0,1,  1,   1]*symScale;
            plot(xVec,yVec,'b');
        end
        if thetaFixed(node)
            xVec = globalPos(node,1) + 0.7*[0,0.5  ,0.866,1]*symScale;
            yVec = globalPos(node,2) + 0.7*[1,0.866,  0.5,0]*symScale;
            plot(xVec,yVec,'g');
        end
        if xMoves(node)
            xVec = globalPos(node,1) + -[0,1,0.6, 0.6,1]*symScale;
            yVec = globalPos(node,2) + -[0,0,0.4,-0.4,0]*symScale;
            plot(xVec,yVec,'r');
        end
        if yMoves(node)
            xVec = globalPos(node,1) + -[0,0,0.4,-0.4,0]*symScale;
            yVec = globalPos(node,2) + -[0,1,0.6, 0.6,1]*symScale;
            plot(xVec,yVec,'b');
        end
        if thetaMoves(node)
            xVec = globalPos(node,1) + 0.7*[0,0.5  ,0.866,1,1.33,1.1518,.665,0]*symScale;
            yVec = globalPos(node,2) + 0.7*[1,0.866,  0.5,0,0,0.665,1.1518,1.3300]*symScale;
            plot(xVec,yVec,'g');
        end
    end
    title('Boundary Conditions')
end %end checkBC

function nodes = rotCenterGroup(limit,rotMag)
%This function adds the nodes in the center groups on the appropriate
%borders to the appropriate fixed boundary condition.
%xFix,yFix, and thetaFix are boolean values for whether those things should
%be fixed.
switch limit
    case 'left'
        centers = centerGroups((abs(centerGroups(:,1)) < tol) ,:);
    case 'right'
        centers = centerGroups(abs(centerGroups(:,1)-macroWidth) <tol,:);
    case 'top'
        centers = centerGroups(abs(centerGroups(:,2)-macroHeight)<tol,:);
    case 'bottom'
        centers = centerGroups(abs(centerGroups(:,2))<tol,:);
end
%Centers now holds the centerGroups that are on the appropriate border.
%Figure out which nodes are in those groups.
temp = centers(:,3:8);
nodes = temp(find(not(centers(:,3:8)==0))); %#ok

%Calculate the displacements that are caused by rotating the center group.
%These are independant from the displacements caused by translation.
rotMat = [cos(rotMag),-sin(rotMag);sin(rotMag),cos(rotMag)];
for i = 1:size(centers,1) %loop over the center groups.
    for j = 3:size(centers,2) %loop over the nodes in each center groups.
        if centers(i,j)
            centerXY = centers(i,1:2);
            xy = globalPos(centers(i,j),:);
            r = (xy-centerXY)';
            uRot(centers(i,j),:) = [(rotMat*r-r)',rotMag];
        end
    end
end
end
%If this is a displacement controlled boundary condition set the free force
%vector.
if not(exist('f_f')==1)
    f_f = zeros(length(free),1);
end

end