##
algb=/home/younguj/repo/abaqusPy/umats/lib/algb.f
lib_write=/home/younguj/repo/abaqusPy/umats/lib/lib_write.f
lib=/home/younguj/repo/abaqusPy/umats/lib/lib.f
is=/home/younguj/repo/abaqusPy/umats/lib/is.f
cnv=/home/younguj/repo/abaqusPy/umats/lib/cnv.f

f2py -c -m yld2000 yld2000_2d.f $algb $lib_write $lib $is $cnv -DF2PY_REPORT_ON_ARRAY_COPY=1.
