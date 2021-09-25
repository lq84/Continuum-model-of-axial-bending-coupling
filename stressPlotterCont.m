function stressPlotterCont(u,IEN,coords,D)
%This function will plot the six stress plots for a continuum problem.
%This will give me a chance to examine what makes the solutions tick.

%Initialize the figure and the the subplots

nn = length(u)/3; ne = size(IEN,1);
%% Initialize variables for recalculating the shape matrices.
%etaRow = [-sqrt(3/5),0,sqrt(3/5)]; etaRow = [etaRow,etaRow,etaRow]; %Gauss points
%xiRow =[-sqrt(3/5),-sqrt(3/5),-sqrt(3/5),0,0,0,sqrt(3/5),sqrt(3/5),sqrt(3/5)]; %Gauss points
%weight = [25,40,25,40,64,40,25,40,25]/81;
%etaRow = [-1,0,1]; etaRow = [etaRow,etaRow,etaRow];
%xiRow =[-1,-1,-1,0,0,0,1,1,1];
dofs = zeros(24,1);
centerCoordsMat = zeros(ne,2);
strainMat = zeros(ne,6);
stressMat = zeros(ne,6);

%% Loop over the elements and calculate the stress at the center.
for el = 1:ne
    nodes = IEN(el,:);
    elCoords = coords(nodes,:);
    centerCoord = mean(elCoords);
    centerCoordsMat(el,:) = centerCoord; 
    %Store the coordinates of the center of the element for later plotting
    dofs(1:3:end) = nodes*3-2;
    dofs(2:3:end) = nodes*3-1;
    dofs(3:3:end) = nodes*3;
    [B, ~] = BMatrixQ(0,0,coords(nodes,:));
    uElement = u(dofs);
    strainMat(el,:) = (B*uElement)';
    stressMat(el,:) = (D*B*uElement)';
end

%% Actually plot that stuff.
figure; 
%X and Y are vectors that show the unique nodal positions.
X = unique(centerCoordsMat(:,1));
Y = unique(centerCoordsMat(:,2));
%Figure out what the contours should be
%lowerBound = min(min(strainMat));
%upperBound = max(max(strainMat));
%cont = linspace(lowerBound,upperBound,14);
strainsTitle = cell(6,1);
strainsTitle{1} = 'epsilon11';
strainsTitle{2} = 'epsilon22';
strainsTitle{3} = 'epsilon12';
strainsTitle{4} = 'epsilon21';
strainsTitle{5} = 'phi3,1';
strainsTitle{6} = 'phi3,2';
for i = 1:6 %Arrange the data so that it can be plotted in a contour.
    subplot(2,3,i);
    
    Z = zeros(length(X),length(Y));
    for j = 1:length(X)
        %Loop over the data and put it into the right spot in the z matrix.
        Z(j,:) = strainMat(X(j)==centerCoordsMat(:,1),i)';
    end
    [~,h] = contourf(X,Y,Z',60);
    set(h,'EdgeColor','none');
    title([strainsTitle{i},' Min=',num2str(min(min(Z))),' Max=',num2str(max(max(Z)))]);
end
%Do it all over again for stress.
figure;
title('Six Stress Components')
%Figure out what the contours should be

strainsTitle = cell(6,1);
strainsTitle{1} = 'sigma11';
strainsTitle{2} = 'sigma22';
strainsTitle{3} = 'sigma21';
strainsTitle{4} = 'sigma12';
strainsTitle{5} = 'm13';
strainsTitle{6} = 'm23';
for i = 1:6 %Arrange the data so that it can be plotted in a contour.
    subplot(2,3,i);
    
    Z = zeros(length(X),length(Y));
    for j = 1:length(X)
        %Loop over the data and put it into the right spot in the z matrix.
        Z(j,:) = stressMat(X(j)==centerCoordsMat(:,1),i)';
    end
    [~,h] = contourf(X,Y,Z',60);
    set(h,'EdgeColor','none');
    title([strainsTitle{i},' Min= ',num2str(min(min(Z))),' Max= ',num2str(max(max(Z)))]);
end