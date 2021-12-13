## Conduct a smoke test for bauschinger module
import os,subprocess
import MP.lib.temp

def gen_fo():
    fn = MP.lib.temp.gen_tempfile()
    return open(fn,'w')


def compile_f2py(files):
    """
    Run f2py -m -c dum
    """
    cmd = 'f2py -m -c dum'.split()
    for f in files:
        cmd.append(f)
    stdout = gen_fo(); stderr = gen_fo()
    p = subprocess.Popen(cmd,stdout=stdout,stderr=stderr)
    p.wait()
    iflag=p.poll()
    stdout.close();stderr.close()

    print('iflag:',iflag)
    print('stdout:',stdout.name)
    print('stderr:',stderr.name)
    return iflag

def main():
    ## Fortran source code Compile test
    print('** Bauschinger/latent/hah compile test')
    files=['bauschinger.f','latent.f','hah.f']
    iflag = compile_f2py(files)
    if iflag==0: print('* success')
    else: print('* fail')

    ## hah.f test
    print('** hah_test.f compile test')
    files=['hah_test.f']
    iflag = compile_f2py(files)
    if iflag==0: print('* success')
    else: print('* fail')

if __name__=='__main__':
    main()
