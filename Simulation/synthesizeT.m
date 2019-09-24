function [TGT,TFov]=synthesizeT(NS,theta)

%SYNTHESIZET   Synthesizes a set of rigid transforms affecting the 
%   different shots.
%   [TGT,TFOV]=SYNTHESIZET(NS,THETA) synthesizes the parameters (rotations 
%   in degrees, translations in pixels generically)
%   * NS is a vector with first row specifying the number of shots (or 
%   of transform states), second row specifying the SENSE factors along
%   each direction and third row specifying the number of within shot
%   subdivisions
%   * THETA is a vector with each column specifying the rotation parameter.
%   row 1 is inter-shot rotation, row 2 is inter-shot translation, row 3 
%   is intra-shot rotation, row 4 is intra-shot translation
%   ** TGT is a cell array of transforms
%   ** TFOV is a cell array of transforms to change the geometry in each
%   repeat (for multiple-orientation data)
%

NS(end+1:5,:)=1;
theta(end+1:4,:)=0;
TFov=cell(1,size(NS,2));
for s=1:size(NS,2)%Shots
    TFov{s}=single(zeros([1 1 1 1 1 1 1 NS(7,s) 6]));
    for r=1:NS(7,s)
        if any(NS(4:5,s)<0)
            TFov{s}(:,:,:,:,:,:,:,r,1:2)=(r-1)*abs(permute(NS(4:5,s),[2 3 4 5 6 7 8 9 1]))/NS(7,s);
        else
            TFov{s}(:,:,:,:,:,:,:,r,4)=(r-1)*pi/NS(7,s);
        end
    end
    for V=1:size(theta,2)%Rotation levels
        NT=NS(1,s)/(prod(abs(NS(4:5,s)))*(NS(6,s)^2));
        
        %Cell array of transforms
        TGT{s}{V}=single(zeros([1 1 1 NS(3,s) NT abs(NS(4:5,s))' NS(7,s) 6]));%6th dimension for transform parameters: 1-3 translations / 4-6 rotations
        TGT{s}{V}=dynInd(TGT{s}{V},4,9,pi*bsxfun(@plus,theta(1,V)*(rand([1 1 1 1 NT abs(NS(4:5,s))' NS(7,s)])-0.5),theta(3,V)*(rand([1 1 1 NS(3,s) NT abs(NS(4:5,s))' NS(7,s)])-0.5))/180);
        TGT{s}{V}=dynInd(TGT{s}{V},1:2,9,bsxfun(@plus,theta(2,V)*(rand([1 1 1 1 NT abs(NS(4:5,s))' NS(7,s) 2])-0.5),theta(4,V)*(rand([1 1 1 NS(3,s) NT abs(NS(4:5,s))' NS(7,s) 2])-0.5)));      
        TGT{s}{V}=bsxfun(@minus,TGT{s}{V},multDimMea(TGT{s}{V},4:8));%Due to considerations above ec (23) in L. Cordero-Grande et al, "Sensitivity encoding for aligned multishot magnetic resonance reconstruction," IEEE TCI, 2(3):266–280.
        NT=size(TGT{s}{V});
        TGT{s}{V}=reshape(TGT{s}{V},[NT(1:3) 1 prod(NT(4:5)) NT(6:9)]);%Grouping of motion states
    end
end
