function [coords,IEN] = altHalfspaceMesh()
%This function creates a finite element mesh for the micropolar and
%classical models. It is much denser near the halfspace load and less dense
%elsewhere.

global nelx nn macroWidth macroHeight
%% Initialize
if mod(nelx,2)
    nelx = nelx+1;
end
nodePerEdge = nelx+1;
%The region is divided into two main parts, the lower corner quad and the
%rest. The lower corner is much more finely meshed than the rest of it. The
%bottom left element has edge lengths that are 1/biasRatio as long as the
%edges of the elements on the outside of the corner quad. Those elements
%are 1/biasRatio as long as the ones on the outside of the whole square.
coords = zeros((nelx+1)*(nelx/2+1)+nelx^2/4,2);
IEN = zeros(nelx^2*3/4,4); % Index of element nodes (IEN)
biasRatio = 5;
cornerFrac = .1;

%Create a vector called biasVec to store a way to space points between two
%nodes in space.
%First store the element lengths then convert them to node positions.
biasVec = linspace(1,biasRatio,nelx/2);
for i = (nelx/2+1):-1:1
    biasVec(i) = sum(biasVec(1:(i-1)));
end
biasVec = biasVec'/biasVec(end);

%The bottom left half of the space can be set up fairly easily. Loop over
%the rows and create a corner quad node set and a connecting outside node
%set.
%Start with a set of reference coordinates for the corner edge, the
%right outer edge and the right edge of the corner quad.
cornerVec = biasVec*cornerFrac; 
outVec = linspace(0,1,nelx/2+1);
%cornerVec and outVec are both strictly y coordinates. The x coordinates
%are either constant zero or constant 1.
quadVecX = linspace(cornerFrac,cornerFrac*1.1*cosd(45),nelx/2+1);
quadVecY = linspace(0,cornerFrac*1.1*sind(45),nelx/2+1);
%quadVecX and Y are the coordinates of the right side of the corner
%quadralateral.
count = 1;
for row = 1:(nelx/2+1)
    %Get the corner quad vectors first.
    biasFrac = (row-1)/(nelx/2);
    tempBiasVec = biasVec*(1-biasFrac)+biasFrac*linspace(0,1,nelx/2+1)';
    tempX = 0+(quadVecX(row)-0)*tempBiasVec;
    tempY = cornerVec(row)+(quadVecY(row)-cornerVec(row))*tempBiasVec;
    %Now do the outside area.
    tempX = [tempX(1:end-1);quadVecX(row)+(1-quadVecX(row))*biasVec];
    tempY = [tempY(1:end-1);quadVecY(row)+(outVec(row)-quadVecY(row))*biasVec];
    coords(count:count+nelx,:) = [tempX,tempY];
    count = count+nodePerEdge;
end
%quadCornerNode stores the node number of the corner of the quad area for
%later.
qCN = count-nelx/2-1;
%plot(coords(:,1),coords(:,2),'x')
%% Set up the IEN
%This next section sets up the IEN for the section that just had its points
%generated. This is copy pasted from the rectangular mesh generator. The
%section generated so far can be transformed into a grid mesh simply by
%moving nodes around.
count = 1; %In this loop count stores the element number.
numNodesInRow = nelx+1;
% Each row, so nely # of row
for row = 1:nelx/2
    rowMultiplier = row-1;
    % Each column, so nelx # of row
    for col= 1:nelx
        IEN(count,:)=[rowMultiplier*    numNodesInRow+col, ...
                      rowMultiplier*    numNodesInRow+col+1, ...
                     (rowMultiplier +1)*numNodesInRow+col+1, ...
                     (rowMultiplier +1)*numNodesInRow+col];
        count = count+1;
    end
end
%BOOM! Dun. With that at least.
%% Now make the rest of the mesh. 
%This is a differently oriented grid. It is
%solely from the top edge of the quad to the top edge of the square.
%swap quadVecX and quadVecY because we are looking at the other edge of the
%quad.
%plot(coords(:,1),coords(:,2),'x')
temp = quadVecX; quadVecX = quadVecY; quadVecY = temp;
%Now reverse the order of both so that we can start our process where we
%left off without reverse counting shenanigans.
quadVecX = quadVecX(end-1:-1:1); quadVecY = quadVecY(end-1:-1:1); 
%Also reverse the order of outVec
outVec = outVec(end-1:-1:1);
countIEN = find(IEN(:,1)==0,1)-1;
count = find(and(coords(:,1)==0,coords(:,2)==0),2);
count = count(2);
for row = 1:(nelx/2)
    %Now do the outside area.
    %plot(coords(:,1),coords(:,2),'x')
    tempX = quadVecX(row)+(outVec(row)-quadVecX(row))*biasVec;
    tempY = quadVecY(row)+(1-quadVecY(row))*biasVec;
    coords(count:(count+nelx/2-1),:) = [tempX(2:end),tempY(2:end)];
    count = count+nelx/2;
    %Lets do the IEN while we are at it.
    
    if row==1
        IEN(countIEN+(1:nelx/2),:) = ...
            [qCN:qCN+nelx/2-1;
            (qCN:qCN+nelx/2-1)+1;
             qCN+nelx/2+1:qCN+nelx;
            [qCN-1,qCN+nelx/2+1:qCN+nelx-1]]';
        countIEN = nelx/2+countIEN;
    else
        temp = [0:nelx/2-1;1:nelx/2;nelx/2+1:nelx;[-1,nelx/2+1:nelx-1]]';
        temp = temp+qCN+(row-1)*nelx/2;
        temp(1,[1,4]) = [qCN-(row-1),qCN-(row)];
        IEN(countIEN+(1:nelx/2),:) = temp;
        countIEN = nelx/2+countIEN;
    end
end
coords(:,1) = coords(:,1)*macroWidth;
coords(:,2) = coords(:,2)*macroHeight;
nn = size(coords,1);
%% Begin visual check
%To use this visual check put a breakpoint between the two plot lines.
%figure
%hold on

%% Convert the mesh from linear to quadratic

IEN = [IEN,zeros(size(IEN,1),size(IEN,2))];
nextNodeNum = size(coords,1);
coords = [coords;zeros(nn,2)];

%Loop over the elements
for el = 1:size(IEN,1);
    %hold on
    %xy = mean(coords(IEN(el,1:4),:));
    %text(xy(1),xy(2),num2str(el))
    for sideNode = 5:8
        %Check to see if there is already a created node in this spot.
        if IEN(el,sideNode)==0 %there is not a created node
            %Create a node
            nextNodeNum = nextNodeNum+1;
            IEN(el,sideNode) = nextNodeNum;
            
            %Search for the other spots that node goes in the IEN
            %Figure out which nodes the new node is between.
            switch sideNode
                case 5
                    betweenNodes = IEN(el,[1,2]);
                case 6
                    betweenNodes = IEN(el,[2,3]);
                case 7
                    betweenNodes = IEN(el,[3,4]);
                case 8
                    betweenNodes = IEN(el,[4,1]);
                otherwise
                    error('What the heck? side node is bad')
            end
            coords(nextNodeNum,:) = mean(coords(betweenNodes,:));
            
            %Find other elements with those two nodes
            %There should only be one besides the current element.
            temp1 = IEN(:,1:4)==betweenNodes(1);
            temp2 = IEN(:,1:4)==betweenNodes(2);
            otherEl = find(sum(temp1+temp2,2)==2);
            otherEl = otherEl(not(otherEl==el));
            if length(otherEl)>1
                error('other element appears to be more than a single');
            end
            if not(isempty(otherEl))
                %Figure out which node this is supposed to be and assign
                %the node to it.
                
                betweenNodes = [find(IEN(otherEl,1:4)==betweenNodes(1)),...
                                find(IEN(otherEl,1:4)==betweenNodes(2))];
                betweenNodes = sort(betweenNodes);
                if isequal(betweenNodes,[1,2])
                    IEN(otherEl,5) = nextNodeNum;
                elseif isequal(betweenNodes,[2,3])
                    IEN(otherEl,6) = nextNodeNum;
                elseif isequal(betweenNodes,[3,4])
                    IEN(otherEl,7) = nextNodeNum;
                elseif isequal(betweenNodes,[1,4])
                    IEN(otherEl,8) = nextNodeNum;
                end
            end
        end
    end
end
%% Visual check of mesh This seems to be buggy
%if not(isempty(strfind(P.plots,'IEN')))
% figure; hold on;
% plot(0,0); xlim([0,max(coords(:,1))]); ylim([0,max(coords(:,2))]);
% nodeOrder = [1,5,2,6,3,7,4,8,1];
% for e = 1:size(IEN,1)
%     nodes = IEN(e,:);
%     center = mean(nodes);
%     nodes = coords(nodes(nodeOrder),:);
%     plot(nodes(:,1),nodes(:,2),'r')
%     for n = 1:8
%         text(nodes(n,1),nodes(n,2),num2str(IEN(e,n)));
%     end
%     plot(nodes(:,1),nodes(:,2),'b')
% end
% 1;
%end