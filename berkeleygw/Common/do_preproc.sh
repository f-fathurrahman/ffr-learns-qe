for fil in *.f90
do
  BASNAM=`basename $fil .f90`
  cpp -C -nostdinc -I.  -DGNU  -DUSEFFTW3 $fil > PROCESSED/$BASNAM.f90
done
