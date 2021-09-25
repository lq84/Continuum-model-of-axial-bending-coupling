function P = defaultRun(varargin)
%This function can be run with no input arguments to run a single triplet.
%The parameters for the triplet are defined within the script. If the
%script is called with any input arguments, it will return a structure P.
%This structure can be input into any of the three FEA solvers and they
%will run that problem. If any function is called without input arguments
%they will call this function and use its output in place of input
%arguments.

global L

%P is an input/output structure. It is capable of storing all of the
%parameters and results for a given problem.
P = struct('numCellsHorz',[],'numCellsVert',[],...%Major problem parameters
           'cellHeight',[],'cellWidth',[],'ligThick',[],...%Major problem parameters
           'latticeType',[],'probType',[], ... %Major problem parameters
           'offset',[],...
           'nelx',[],'nely',[],'plots',[], ...
           'AR',[],... %end of solver parameters
           'strainEnergyCont',[], ... %Raw micropolar results
           'strainEnergyClass',[],... 'UCStrainEnergyClass',[], ...
           'globalStrainEnergyError',[],...
           'globalClassEnergyError',[],...
           'ligLength',[],... %Calculated parameters that describe the lattice.
           'useNewProps',true,... %If true, the continuum solver will use the new correct properties from generalized continuum modeling.
           'MPEffectSize',[]);
P.numCellsHorz = 4; %Number of cells in the horizontal direction.
P.numCellsVert = 4; %Number of cells in the vertical direction.
P.latticeType = 'diamondOld'; %This sets the type of topology
%Common options include:
%'triangle','square','hexagon','mixedTriAold','mixedTriBold','mixedTriAnew','mixedTriBnew'
P.probType = 'stretchY'; %This sets the category of boundary conditions.
%Common options include:
%'stretchX','stretchY','stretchYfree','transverseY','transverseX','transverseYhinge','curveX','curveY','halfspace',

%% Figure out how to make the smaller macro-dimension 100
%The following few lines of code set the cellWidth and cellHeight to make
%the smaller macro-dimension to be 100.
[cellAR,thicknessRatio] = isotropic(P.latticeType,1);
macroAR = cellAR*P.numCellsVert/P.numCellsHorz;
if macroAR<=1
    P.cellHeight = L/P.numCellsVert;
    P.cellWidth = P.cellHeight/cellAR;
else
    P.cellWidth = L/P.numCellsHorz;
    P.cellHeight = P.cellWidth*cellAR;
end
%The ligament thickness is going to be 1 tenth of the length.
P.ligThick = P.cellWidth*thicknessRatio;
P.offset = [0,0]; %offset is given in fractions of a unit cell.
%P.offset = [0,1/3]; %Default offset for hexagon
%P.posDef = true; 
%Kumar's paper 'Generalized continuum modeling...' presented two sets of
%material properties. One with gamma>0, one with gamma<0. If P.posDef=true,
%the formulas for gamma>0 are used. This functionality was not discussed in
%the paper.
P.plots = 'displacement';
%P.plots = 'stress,displacement,IEN,blackLattice';
%This line will dictate which plotting options are used as the code runs.
%To include more than one plot option, have a single string where the
%options are separated by commas.
%Options: 'stress,displacement,IEN,blackLattice'
%WARNING: If Matlab opens up hundreds of plots, you can bring your computer
%to a grinding halt.

P.nelx = 8;%NUmber of finite elements in the x and y directions. Default is 35.
P.nely = 8;

%If defaultRun has no input arguments, run a trio of simulations.
if isempty(varargin)
    P = latticeSolver(P);
    P = continuumSolverQ(P);
    P = classicalSolverQ(P);
    P = myError(P);
    disp(['            MP global strain energy error = ',num2str(P.globalStrainEnergyError*100),'%'])
    disp(['         Class global strain energy error = ',num2str(P.globalClassEnergyError*100) ,'%'])
    disp(['MP Effect size global strain energy error = ',num2str(P.MPEffectSize*100) ,'%'])
end

end