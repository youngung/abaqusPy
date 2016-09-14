"""
Compiling test for umats.
"""


## Umat header that includes program
UmatHeader="""
program main
write(*,*)'Hi, someone testing a umat'
end
"""

UmatSubrXit="""
subroutine xit
write(*,*) 'Hello, Subroutine xit was called.'
return
end subroutine xit
"""


def string2IndentedLines(string='',nind=6):
    """
    Get string and split them into lines contained in a list.
    For each line, add indents specified by nindent.
    The new lines are returned as a list.

    Arguments
    ---------
    string
    nind = 6
    """
    lines=string.split('\n')
    indent=' '*nind ## spacing
    for i in xrange(len(lines)):
        lines[i]='%s%s'%(indent,lines[i])
    return lines

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

    nind=6 ## number of indentation used in the fortran fixed form...
    mxLineNo=20 ## max line number to comment out 'include' command

    ## UMAT subroutine Lines in list
    umatSubrLines = open(umat,'r').read().split('\n') ## list
    ## - find <include 'aba_param.inc'> and comment that line out.

    LinesUnderExamine = umatSubrLines[:mxLineNo]
    ## detect aba_param.inc and comment that line out.
    for i in xrange(len(LinesUnderExamine)):
        line = LinesUnderExamine[i]
        if line[nind:nind+7].upper()=='INCLUDE':
            LinesUnderExamine[i]='c %s'%LinesUnderExamine[i]
            print 'Found INCLUDE line and commented out'
        else: pass


    ## UMAT Head Lines in list
    umatHeadLines=string2IndentedLines(string=UmatHeader,nind=nind)

    ## UMAT xit Lines in list
    umatXitLines=string2IndentedLines(string=UmatSubrXit,nind=nind)

    umatProgram = umatHeadLines+umatSubrLines+umatXitLines

    s=''
    for i in xrange(len(umatProgram)):
        s='%s%s\n'%(s,umatProgram[i])

    ## construct a temp file...
    import MP.lib.temp
    # fn = MP.lib.temp.gen_tempfile(prefix='UMAT-TEST',ext='f',tmp='/tmp/younguj')
    fn = MP.lib.temp.gen_tempfile(prefix='UMAT-TEST',ext='f',tmp='/scratch1/younguj/')
    print 'fn:',fn
    with open(fn,'w') as fo:
        fo.write(s)

    if verbose:
        print s

    ## compiling option
    cmd = 'gfotran %s'%(fn)

## Command line usage will be much easier to use.
if __name__=='__main__':
    import argparse, os

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-umat',type=str,help='UMAT file (full-pathed)',
        default='/home/younguj/repo/abaqusPy/umats/el/iso.f')

    args = parser.parse_args()
    main(umat=args.umat)
