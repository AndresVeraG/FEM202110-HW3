**discussion for homework 3.3**

First of all, note that the model has lots of flaws. The most important is that the model does not recognize or distinguish between the estremes and the plates, so the stresses and displacements are asociated with the right extreme of the model. Aditionally, it can be noted that the stress averaging done by the quad4 elements done with tributary areas are not so easiy implemented in a quad 9 element , so that a non invertible matrix in the firts quadrant unabled to continue the analisis. The results obtained of this model are not complete and have some issues to be bettered but can be noted from hipothesis of the knowledge of FEM that:
1. A good implementation of QUAD9 elements can distribute better the stressses and displacements with less interference between the elements.
2. a good implemnetration of stress averaging can reduce the overlock of the nodes show better representation of the stresses in the elemnets because the forces are more uniformelly distributed in the lements. 
3. In the cases studied, the case might be that the shape of the quad 9, show the hipothetical cases described before.

Bellow it can be seen the images represented in the cases shown above. 


