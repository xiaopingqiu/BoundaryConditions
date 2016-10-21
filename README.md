In the context of finite volume discretization, the implemantation of no slip wall boundary condition should guarantee that the wall shear stress equals to viscosity times the normal gradient of velocity parallel to the wall. But In OpenFOAM, no slip wall condition is simplified to a Dirichlet boundary condition, which may introduce important errors.


This repository implemented the "real no slip wall boundary condition" for OpenFOAM-2.3.x and OpenFOAM-4.x. Please switch to the appropriate branch according to the version of your OpenFOAM. 


Warning: the codes are not intensively tested, use them at your own risks. If you find bugs, please leave a message. 


**Reference**:

[The Finite Volume Method in Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM® and Matlab®](http://www.springer.com/us/book/9783319168739)
