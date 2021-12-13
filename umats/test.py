"""
Compiling test for umats.
"""

def loadModules():
    """
    Load modules necessary in Palmetto
    """
    import os
    modules=['abaqus/6.14','intel/16.0']

    cmd='module load'
    for module in modules:
        cmd ='%s %s'%(cmd,module)

    print('cmd:')
    print(cmd)
    os.system(cmd)


def main(umat='/home/younguj/repo/abaqusPy/umats/el/iso.f',
         verbose=True):
    """
    Generate a temporary file that wraps the given <umat>
    to test 'debug' the file using external fortran compiler
    to minimize the syntax errors...

    Arguments
    ---------
    umat='/home/younguj/repo/abaqusPy/umats/el/iso.f'
    verbose=True
    """
    import subprocess

    loadModules()

    print('-'*40)
    print('UMAT file name')
    print(umat)
    print('-'*40)

    nind=6 ## number of indentation used in the fortran fixed form...
    mxLineNo=20 ## max line number to comment out 'include' command

    ## UMAT subroutine Lines in list
    umatSubrLines = open(umat,'r').read().split('\n') ## list
    ## Find <include 'aba_param.inc'> and comment that line out.
    LinesUnderExamine = umatSubrLines[:mxLineNo]
    ## Detect aba_param.inc and comment that line out.
    for i in range(len(LinesUnderExamine)):
        line = LinesUnderExamine[i]
        if line[nind:nind+7].upper()=='INCLUDE':
            umatSubrLines[i]='C %s'%LinesUnderExamine[i]
        else: pass

    ## path to fake main
    fnFakeMain = 'fakemain.f'
    mainLines = open(fnFakeMain, 'r').read().split('\n')

    umatProgram = mainLines+umatSubrLines

    s=''
    for i in range(len(umatProgram)):
        s='%s%s\n'%(s,umatProgram[i])

    ## construct a temp file...
    import MP.lib.temp
    fnTempUmat = MP.lib.temp.gen_tempfile(prefix='UMAT-TEST',ext='f',tmp='/scratch1/younguj/')
    print('fn:',fnTempUmat)
    with open(fnTempUmat,'w') as fo:
        fo.write(s)
        pass

    if verbose: print(s)

    ## compiling option
    # compiler='gfortran'  ## GNU fortran compiler
    compiler='ifort'     ## intel fortran compiler

    cmd='%s %s > compile.log'%(compiler,fnTempUmat)
    print('cmd:')
    print(cmd)
    iflag=os.system(cmd)
    if iflag!=0:
        print('Compile failed')
    else:
        cmd='./a.out'
        iflag_run=os.system(cmd)
        if iflag_run!=0:
            print('Running error')
        else:
            print('Running success')
            print('printout fort.7 the msg file')
            cmd='cat fort.7'
            os.system(cmd)

            ## add umat
            path_home = os.getcwd()
            path_run = os.path.split(umat)[0]
            os.chdir(path_run)
            cmd ='abaqus make library=%s'%(umat)
            os.system(cmd)
            os.chdir(path_home)


    # cmd=[compiler,fnTempUmat]
    # subprocess.check_output(cmd)

## Command line usage will be much easier to use.
if __name__=='__main__':


    ## add more libraries as progress
    libraries=['/home/younguj/repo/abaqusPy/umats/el/iso.f',
               '/home/younguj/repo/abaqusPy/umats/epl/mises.f']
    print('-'*72)
    print('Available libraries')
    for lib in libraries:
        print(lib)
    print('-'*72)



    ## Paths to UMAT files
    default_UMAT_FN='/home/younguj/repo/abaqusPy/umats/el/iso.f'
    # default_UMAT_FN='/home/younguj/repo/abaqusPy/umats/epl/mises.f'

    import argparse, os

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-umat',type=str,help='UMAT file (full-pathed)',
        default=default_UMAT_FN)
    parser.add_argument('-v',action='store_true',
                        help='Verbose flag')

    args = parser.parse_args()
    main(umat=args.umat,verbose=args.v)
