## Running one.py without gui


## uniaxial tension along global X at 3 different directions
abaqus cae noGUI=one.py -- -iwait -iload 0 -nth 3 -umat /home/younguj/repo/abaqusPy/umats/epl/epl.f
## uniaxial tension along global Y at 3 different directions
abaqus cae noGUI=one.py -- -iwait -iload 3 -nth 3 -umat /home/younguj/repo/abaqusPy/umats/epl/epl.f
abaqus cae noGUI=one.py -- -iwait -iload 2 -nth 3 -umat /home/younguj/repo/abaqusPy/umats/epl/epl.f
