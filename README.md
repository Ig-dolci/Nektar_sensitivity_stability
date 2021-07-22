
Nektar ++ is an open-source software framework designed to support the
development of high-performance scalable solvers for partial differential
equations (PDEs) using the spectral/hp element method.

The original software and User Guide are available for download from
<http://www.nektar.info/>.



Tutorials
---------
A number of tutorials are available, designed to walk the user through the
basics of spectral/hp element methods, through the use of individual solvers and
performing specific types of calculations.

The tutorials are available from <http://doc.nektar.info/tutorials/latest>.

Linear stability analysis for flow around an elastically-mounted rigid structure free to translate.
----------
This repository is focused on adapting the Nektar++ to provide computations of sensitivity and stability analyses for fluid flow system, and fluid-structure interaction system. In the last case, the structure is rigid and only able to be translated.

At the demo repository, there is an example for the user to run a case of linear stability analysis for a flow around an elastically-mounted cylinder free to oscillate on the cross flow direction (y-axis). To do so, follow the steps presented below.

Step 1 - Run basefield.xml to obtain the steady base. In the current example, Reynolds number 33 is set.

Step 2 - As introduced recently in the paper https://doi.org/10.1017/jfm.2020.685, it is necessary to get the gradient of the base flow aerodynamic forces at the wall of the structure. To do it, run the command 
FieldConvert -m gradient basefield.xml basefield.fld basefield_grad.fld, next the command FieldConvert -m extract:bnd=1 basefield.xml basefield_grad.fld boundary_extration.fld

Step 3 - Taking the fields of basefield.fld and boundary_extration_b1.fld, the last step is to run the direct.xml to obtain the modes and their eigenvalues computed by ARPACK library. In the current example, the linear analysis is carried out for mass ratio 50,  damping parameter equal zero and reduced velocity 8. Definitions of these parameters can be found at the paper https://doi.org/10.1017/jfm.2019.754.
