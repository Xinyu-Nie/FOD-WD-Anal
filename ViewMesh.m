function p = ViewMesh(ptVec,trgVec,color)
%p = ViewMesh(ptVec,trgVec,color)
%
%Plot a mesh 
%
%Yonggang Shi
%May 2015

if nargin==2
%     color =[18 219 200]./255; 
     color = [1.0 0.78 0.67];
end

%figure
p = patch('Vertices',ptVec, 'Faces',trgVec);
set(p,'FaceColor',color);
set(p,'edgecolor','none');
daspect([1 1 1])
view(3)
axis tight
 camlight
lighting gouraud
p.AmbientStrength = 0.5;
p.DiffuseStrength = 0.5;
%lightangle(-45,30)