AbaqusPy
--------
Collection of Python scripts/UMATS used for my research activities based on Abaqus.
This package is under development, will add more material constitutive models.


Libraries
---------
* return_mapping
  forward euler return-mapping subroutine
* yld2000-2d
* hill48


Python Examples
---------------
* One element example

 Running one element test is as easy as running below in the command line:

$> abaqus cae noGUI=one.py


The above will use your UMAT (e.g., epl.f as included in this package) and run for
0.05 strain for 0,45, and 90 degrees from RD.

The above will generate a number of "*.odb" files.
Below command is used to extract some useful data to text files

$> abaqus cae noGUI=onePP.py

which will read any '*.odb' file in the cwd. Then it will extract UVARM and save some data to
ASCII file which can be visualized by <oneElementPlotter.ipynb>


* Full size E8 element tensile tests.
 Multi element tests using a full-size E8 standard uniaxial tensile bar.

 $> abaqus cae noGUI=uni.py

* 2D draw-bending / springback simulation using epl.f - with yld2000-2d

 $> abaqus cae noGUI=spring.py

 This will fetch springback_std_both.inp and copy it and modify it to use
 yld2000-2d as the constitutive model.


User Materials
--------------
epl.f : UMAT based on elasto-plastic constitutitve model

epl.f runs currently for shell element and accepts material orientations.
The standard field varibles will be written in reference to each
material orientation pertaining to element. Field variables in reference to
global coordinates are calculated in UVARM, which gives

UVARM1 : ee11
UVARM2 : ee22
UVARM3 : ee12

UVARM4 : pe11
UVARM5 : pe22
UVARM6 : pe12

UVARM7 : s11
UVARM8 : s22
UVARM9 : s12

Youngung Jeong @ Changwon National University
School of Materials Science and Engineering
yjeong@changwon.ac.kr
