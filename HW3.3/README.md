**discussion for homework 3.3**

First of all, note that the model has lots of flaws. The most important is that the model does not recognize or distinguish between the estremes and the plates, so the stresses and displacements are asociated with the right extreme of the model. Aditionally, it can be noted that the stress averaging done by the quad4 elements done with tributary areas are not so easiy implemented in a quad 9 element , so that a non invertible matrix in the firts quadrant unabled to continue the analisis. The results obtained of this model are not complete and have some issues to be bettered but can be noted from hipothesis of the knowledge of FEM that:
1. A good implementation of QUAD9 elements can distribute better the stressses and displacements with less interference between the elements.
2. a good implemnetration of stress averaging can reduce the overlock of the nodes show better representation of the stresses in the elemnets because the forces are more uniformelly distributed in the lements. 
3. In the cases studied, the case might be that the shape of the quad 9, show the hipothetical cases described before.

Bellow it can be seen the images represented in the cases shown above. 

**coarse**

![Grueso](https://user-images.githubusercontent.com/69157203/120057965-71bd2700-c015-11eb-9a10-3dc53c9045ae.png)

![sig x - grueso](https://user-images.githubusercontent.com/69157203/120057968-77b30800-c015-11eb-9563-bc61dbce3b37.png)

![desp x -grueso](https://user-images.githubusercontent.com/69157203/120057971-7a156200-c015-11eb-84cc-f56c7da699f7.png)


**medium** 

![medio](https://user-images.githubusercontent.com/69157203/120057978-7da8e900-c015-11eb-964c-48b627dc4a38.png)
![sig x - med](https://user-images.githubusercontent.com/69157203/120057986-81d50680-c015-11eb-962c-b6f81a0fe92d.png)
![desp x - med](https://user-images.githubusercontent.com/69157203/120057988-839eca00-c015-11eb-9ccc-ac8d16cdbd6b.png)


**fine** 

![Figure 2021-05-29 001526](https://user-images.githubusercontent.com/69157203/120057992-8994ab00-c015-11eb-86fc-ba962a8b5d3a.png)
![sig x - fine](https://user-images.githubusercontent.com/69157203/120057993-8bf70500-c015-11eb-8006-05be810fb418.png)
![desp x - fine](https://user-images.githubusercontent.com/69157203/120057995-8dc0c880-c015-11eb-8f76-16193571ea26.png)


as can be seen, the refinement in the distribution of stresses and displacements are more unifromelly and less disperse as the mesh refines. As can be noted in HW3., the differences are quite beeg in meshs. So It remains pending the good modelling to det the differnet distributions of stresses and displacemnets in the rest of the object.
