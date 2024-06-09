function Weight =  LinearWeights(point)

T=point;
Weight(:,8)=T(:,1).*T(:,2).*T(:,3);
Weight(:,7)=T(:,1).*T(:,2).*(1-T(:,3));
Weight(:,6)=T(:,1).*(1-T(:,2)).*T(:,3);
Weight(:,5)=T(:,1).*(1-T(:,2)).*(1-T(:,3));
Weight(:,4)=(1-T(:,1)).*T(:,2).*T(:,3);
Weight(:,3)=(1-T(:,1)).*T(:,2).*(1-T(:,3));
Weight(:,2)=(1-T(:,1)).*(1-T(:,2)).*T(:,3);
Weight(:,1)=(1-T(:,1)).*(1-T(:,2)).*(1-T(:,3));
end