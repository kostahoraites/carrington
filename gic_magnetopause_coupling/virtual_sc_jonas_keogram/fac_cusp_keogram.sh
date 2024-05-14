var='proton/vg_rho'  #'vg_j_parallel'
startfile=501
endfile=510 # 1598
bulkpath="/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1"
bulkprefix="bulk1"
outputdir="./"
outputname="test.png"
intpol=1
filt=0
#op='x'
cmap='bwr'
pointfile='/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/virtual_sc_jonas_keogram/fg_innerboundary_btrace_FHAFGB.txt'

#    jplots(
#        var=args.var,
#        fnr0=args.fnr[0],
#        fnr1=args.fnr[1],
#        start_coords=args.startpoint,
#        end_coords=args.endpoint,
#        dr=args.dr,
#        bulkpath=args.bulkpath,
#        bulkprefix=args.bulkprefix,
#        outputdir=args.outputdir,
#        outputname=args.outputname,
#        intpol=args.intpol,
#        filt=args.filt,
#        op=args.op,
#        cmap=args.cmap,
#        pointfile=args.pointfile,
#    )


#python /home/horakons/proj/analysator/scripts/cutthrough_timeseries.py -var $var -fnr 1000 1200 -bulkpath $bulkpath -bulkprefix $bulkprefix -outputdir $outputdir -outputname $outputname -intpol $intpol -filt $filt -op $op -cmap $cmap -pointfile $pointfile
python /home/horakons/proj/analysator/scripts/cutthrough_timeseries.py -var $var -fnr $startfile $endfile -bulkpath $bulkpath -bulkprefix $bulkprefix -outputdir $outputdir -outputname $outputname -intpol $intpol -filt $filt -cmap $cmap -pointfile $pointfile

