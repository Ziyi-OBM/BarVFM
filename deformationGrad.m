function [deformGrad,strainField,XYZcg]=deformationGrad(Node,Displacement,Elements,isVirtualStrain)




%Find strain field
numEle=size(Elements,1);
strainField=zeros(3,3,numEle);
deformGrad=zeros(3,3,numEle);
derivCoeffR=[-1 1 1 -1 -1 1 1 -1];
derivCoeffS=[-1 -1 1 1 -1 -1 1 1];
derivCoeffT=[-1 -1 -1 -1 1 1 1 1];
XYZcg=zeros(numEle,3);

for i=1:numEle
   elemNode=Node(Elements(i,:),:);
   elemDisp=Displacement(Elements(i,:),:);
   
   XYZcg(i,1)=1/8*sum(elemNode(:,1));
   XYZcg(i,2)=1/8*sum(elemNode(:,2));
   XYZcg(i,3)=1/8*sum(elemNode(:,3));
         
    ur=1/8*derivCoeffR*elemDisp(:,1);
    us=1/8*derivCoeffS*elemDisp(:,1);         
    ut=1/8*derivCoeffT*elemDisp(:,1);       
             
    vr=1/8*derivCoeffR*elemDisp(:,2);
    vs=1/8*derivCoeffS*elemDisp(:,2);         
    vt=1/8*derivCoeffT*elemDisp(:,2);  
    
    wr=1/8*derivCoeffR*elemDisp(:,3);
    ws=1/8*derivCoeffS*elemDisp(:,3);         
    wt=1/8*derivCoeffT*elemDisp(:,3);  
    
    xr=1/8*derivCoeffR*elemNode(:,1);
    xs=1/8*derivCoeffS*elemNode(:,1);         
    xt=1/8*derivCoeffT*elemNode(:,1);       
             
    yr=1/8*derivCoeffR*elemNode(:,2);
    ys=1/8*derivCoeffS*elemNode(:,2);         
    yt=1/8*derivCoeffT*elemNode(:,2);  
    
    zr=1/8*derivCoeffR*elemNode(:,3);
    zs=1/8*derivCoeffS*elemNode(:,3);         
    zt=1/8*derivCoeffT*elemNode(:,3);  
    
   pDerivDisp=[ur us ut; vr vs vt; wr ws wt];
   pDerivNode=[xr xs xt; yr ys yt; zr zs zt];            

   eleU=pDerivDisp/pDerivNode;

   if isVirtualStrain==false
   eleF=eleU+eye(3);  
   deformGrad(:,:,i)=eleF; 
   strainField(:,:,i)=1/2.*(eleF'*eleF-eye(3));
   else
       deformGrad(:,:,i)=eleU; 
       strainField(:,:,i)=1/2.*(eleU'+eleU);
   end
   
   
end


end