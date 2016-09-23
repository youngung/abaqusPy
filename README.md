AbaqusPy
--------
Collection of Python scripts/UMATS used for my research activities based on Abaqus.
This package is under development, will add more material constitutive models.


Libraries
---------


Python Examples
---------------
* One element example

 Running one element test is as easy as running below in the command line:

$> abaqus cae noGUI=one.py -- -umat <your user mat file name> -iumat -mxe 0.05 -nth 3

The above will use your UMAT (e.g., epl.f as included in this package) and run for
0.05 strain for 0,45, and 90 degrees from RD.


The above will generate a number of "*.odb" files.
Below command is used to extract some useful data to text files

$> abaqus cae noGUI=onePP.py


* Full size E8 element tensile tests.
 Multi element tests using a full-size E8 standard uniaxial tensile bar.

User Materials
--------------
epl.f : UMAT based on elasto-plastic constitutitve model


Youngung Jeong
youngung.jeong@gmail.com
