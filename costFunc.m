function [ outPut ] = costFunc( inputParam )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% load('TestCaseNeo.mat');
% load('TestCaseNeo_1.mat');
% load('TestCaseNeo_2.mat');
% load('TestCaseNeo_3.mat');
% load('TestCaseNeo_4.mat');
load('TestCaseNeo_5.mat');

% patch('Faces',BoundaryFaces,'Vertices',Node,'FaceAlpha',0)
% axis equal

C1=inputParam(1);
K=inputParam(2);

% patch('Faces',BoundaryFaces,'Vertices',Node,'FaceAlpha',0)
% axis equal

numNode=size(Node,1);
%Find strain field
numEle=size(Elements,1);
strainField=zeros(3,3,numEle);
deformGrad=zeros(3,3,numEle);
for i=1:numEle
    elemNode=Node(Elements(i,:),:);
    elemNodeNew=elemNode+Displacement(Elements(i,:),:);
    eleF1=[     elemNodeNew(2,1)-elemNodeNew(1,1), elemNodeNew(4,1)-elemNodeNew(1,1), elemNodeNew(5,1)-elemNodeNew(1,1);
                elemNodeNew(2,2)-elemNodeNew(1,2), elemNodeNew(4,2)-elemNodeNew(1,2), elemNodeNew(5,2)-elemNodeNew(1,2);
                elemNodeNew(2,3)-elemNodeNew(1,3), elemNodeNew(4,3)-elemNodeNew(1,3), elemNodeNew(5,3)-elemNodeNew(1,3)]'./...
                [elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3);
                 elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3);
                 elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3)];
             
    eleF2= [    elemNodeNew(7,1)-elemNodeNew(8,1), elemNodeNew(7,1)-elemNodeNew(6,1), elemNodeNew(7,1)-elemNodeNew(3,1);
                elemNodeNew(7,2)-elemNodeNew(8,2), elemNodeNew(7,2)-elemNodeNew(6,2), elemNodeNew(7,2)-elemNodeNew(3,2);
                elemNodeNew(7,3)-elemNodeNew(8,3), elemNodeNew(7,3)-elemNodeNew(6,3), elemNodeNew(7,3)-elemNodeNew(3,3)]'./...
                [elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3);
                 elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3);
                 elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3)];   

   eleF=(eleF1+eleF2)./2;

   deformGrad(:,:,i)=eleF; 
   strainField(:,:,i)=1/2.*(eleF'*eleF-eye(3));
end
%

%Find the stress tensor

stressTensor=zeros(3,3,numEle);
strnEngDriv=zeros(3,3,numEle);
for i=1:numEle
    Ii=trace(deformGrad(:,:,i)'*deformGrad(:,:,i));
    lambda1=deformGrad(1,1,i);
    lambda2=deformGrad(2,2,i);
    lambda3=deformGrad(3,3,i);
    jcbian=det(deformGrad(:,:,i));
    
    sigm1=1/jcbian*(2*C1*(lambda1^2-1)+K*log(jcbian));
    sigm2=1/jcbian*(2*C1*(lambda2^2-1)+K*log(jcbian));
    sigm3=1/jcbian*(2*C1*(lambda3^2-1)+K*log(jcbian));
    strnEngDriv(:,:,i)=[sigm1,0,0;0,sigm2,0;0,0,sigm3];
    %stressTensor(:,:,i)=1/jcbian*strnEngDriv(:,:,i)*deformGrad(:,:,i)';
    stressTensor(:,:,i)=strnEngDriv(:,:,i); %%% You double multiplied by F and 1/J
end


%Analyze Virtual Disp 1  u1: x 0 0
u1=[1 0 0; 0 0 0; 0 0 0];

%Find RHS
logicState=Node(BoundaryFaces(:,1),1)==5&...    %Logic array representing
    Node(BoundaryFaces(:,2),1)==5&...           %the faces
    Node(BoundaryFaces(:,3),1)==5&...
    Node(BoundaryFaces(:,4),1)==5;

xFaces=BoundaryFaces(logicState,:);  %Faces on the boundary

num_face=size(xFaces,1);
virDisVec=zeros(4,3,num_face);
virWorkEle=zeros(1,num_face);
for i=1:num_face
    currface=xFaces(i,:);
    virDisVec(:,:,i)=Node(currface,:)*u1;
    virWorkEle(:,i)=mean(virDisVec(:,:,i),1)... %Find the mean of the displacement vector among the four nodes
        *[1;0;0]*EndPressure*(0.2)^2; %dot the pressure,normalize by the face area,times 2 for the other face
end

RHS1= sum(virWorkEle);

%LHS
%Find virtual strain
virStrainField=zeros(3,3,numEle);
virDisplacement=Node*u1;
for i=1:numEle
    elemNode=Node(Elements(i,:),:);
    %elemNodeNew=elemNode+virDisplacement(Elements(i,:),:);
    elemNodeNew=virDisplacement(Elements(i,:),:); %%%
    
    %%% To make the VF method work, they define strain differently. It is
    %%% based on the gradient of the displacement. %I changed the
    %%% definition of elemNodeNew to just be the Virtual Displacement above


     eleF1=[     elemNodeNew(2,1)-elemNodeNew(1,1), elemNodeNew(4,1)-elemNodeNew(1,1), elemNodeNew(5,1)-elemNodeNew(1,1);
                elemNodeNew(2,2)-elemNodeNew(1,2), elemNodeNew(4,2)-elemNodeNew(1,2), elemNodeNew(5,2)-elemNodeNew(1,2);
                elemNodeNew(2,3)-elemNodeNew(1,3), elemNodeNew(4,3)-elemNodeNew(1,3), elemNodeNew(5,3)-elemNodeNew(1,3)]'./...
                [elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3);
                 elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3);
                 elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3)];
             
    eleF2= [    elemNodeNew(7,1)-elemNodeNew(8,1), elemNodeNew(7,1)-elemNodeNew(6,1), elemNodeNew(7,1)-elemNodeNew(3,1);
                elemNodeNew(7,2)-elemNodeNew(8,2), elemNodeNew(7,2)-elemNodeNew(6,2), elemNodeNew(7,2)-elemNodeNew(3,2);
                elemNodeNew(7,3)-elemNodeNew(8,3), elemNodeNew(7,3)-elemNodeNew(6,3), elemNodeNew(7,3)-elemNodeNew(3,3)]'./...
                [elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3);
                 elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3);
                 elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3)];   


             
             eleF=(eleF1+eleF2)./2;
   
   %virStrainField(:,:,i)=1/2.*(eleF'*eleF-eye(3));
   virStrainField(:,:,i)=1/2.*(eleF'+eleF); % This is the defintion they use.
end

%Calculate intergral
internal=zeros(1,numEle);
for i=1:numEle
    internal(i)=sum(sum(stressTensor(:,:,i).*virStrainField(:,:,i)'*(0.2)^3));   
end
LHS1=sum(internal);



%Analyze Virtual Disp 2    u2: 0 1/2y 1/2z
u2=[0 0 0 ; 0 -0.5 0; 0 0 -0.5];


logicState=Node(BoundaryFaces(:,1),2)==1|...    %Logic array representing
    Node(BoundaryFaces(:,2),2)==0|...           %the faces
    Node(BoundaryFaces(:,3),3)==1|...
    Node(BoundaryFaces(:,4),3)==0;

%Find RHS
xFaces=BoundaryFaces(logicState,:);  %Faces on the boundary

num_face=size(xFaces,1);
virDisVec=zeros(4,3,num_face);
virWorkEle=zeros(1,num_face);
for i=1:num_face
    currface=xFaces(i,:);
    virDisVec(:,:,i)=Node(currface,:)*u2;
    virWorkEle(:,i)=mean(virDisVec(:,:,i),1)... %Find the mean of the displacement vector among the four nodes
        *[1;0;0]*EndPressure*(0.2)^2; %dot the pressure,normalize by the face area,times 2 for the other face
end

RHS2= sum(virWorkEle);

%LHS
%Find virtual strain
virStrainField=zeros(3,3,numEle);
virDisplacement=Displacement*u2;
for i=1:numEle
    elemNode=Node(Elements(i,:),:);
    elemNodeNew=virDisplacement(Elements(i,:),:); %% Same changes as above
    eleF1=[     elemNodeNew(2,1)-elemNodeNew(1,1), elemNodeNew(4,1)-elemNodeNew(1,1), elemNodeNew(5,1)-elemNodeNew(1,1);
                elemNodeNew(2,2)-elemNodeNew(1,2), elemNodeNew(4,2)-elemNodeNew(1,2), elemNodeNew(5,2)-elemNodeNew(1,2);
                elemNodeNew(2,3)-elemNodeNew(1,3), elemNodeNew(4,3)-elemNodeNew(1,3), elemNodeNew(5,3)-elemNodeNew(1,3)]'./...
                [elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3);
                 elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3);
                 elemNode(2,1)-elemNode(1,1),elemNode(4,2)-elemNode(1,2),elemNode(5,3)-elemNode(1,3)]';
             
    eleF2= [    elemNodeNew(7,1)-elemNodeNew(8,1), elemNodeNew(7,1)-elemNodeNew(6,1), elemNodeNew(7,1)-elemNodeNew(3,1);
                elemNodeNew(7,2)-elemNodeNew(8,2), elemNodeNew(7,2)-elemNodeNew(6,2), elemNodeNew(7,2)-elemNodeNew(3,2);
                elemNodeNew(7,3)-elemNodeNew(8,3), elemNodeNew(7,3)-elemNodeNew(6,3), elemNodeNew(7,3)-elemNodeNew(3,3)]'./...
                [elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3);
                 elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3);
                 elemNode(7,1)-elemNode(8,1),elemNode(7,2)-elemNode(6,2),elemNode(7,3)-elemNode(3,3)]';   

   eleF=(eleF1+eleF2)./2;

   virStrainField(:,:,i)=1/2.*(eleF'+eleF);
end


%Calculate intergral
internal=zeros(1,numNode);
for i=1:numEle
    internal(i)=sum(sum(stressTensor(:,:,i).*virStrainField(:,:,i)'*0.2^3));   
end
LHS2=sum(internal);


outPut=(LHS1-RHS1)^2+(LHS2-RHS2)^2;




end

