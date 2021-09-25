function [height,thickness] = isotropic(type,width)
    switch type
        case {'triangle','hexagon'}
            height = width*sqrt(3);
        case 'hexachiral'
            height = width/sqrt(3);
        case {'square','mixedTriAold','mixedTriBold','mixedTriAnew','mixedTriBnew','hexSquare','diamondNew','diamondOld','tigerPaw'}
            height = width;
        otherwise
            warning('defaultRun->isometric does not recognize lattice type')
            height = width;
    end
    
    %Figure out ligament length and therefore thickness.
    switch type
        case {'triangle','square'}
            L = width;
        case {'mixedTriAold','diamondOld'}
            L1 = width; L2 = sqrt(2)/2*width;
            L = L1/2+L2/2;
        case 'mixedTriBold'
            L1 = width/2; L2 = width/2*sqrt(2);
            L = L1/2+L2/2;
        case {'mixedTriAnew','diamondNew'}
            L1 = width; L2 = width*sqrt(2);
            L = L1/2+L2/2;
        case 'mixedTriBnew'
            L1 = width; L2 = width/2*sqrt(2);
            L = L1/2+L2/2;
        case {'hexagon'}
            L = width/2/cosd(30);
    end
    thickness = L/10;
end