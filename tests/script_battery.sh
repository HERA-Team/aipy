#! /bin/bash

# Exercise lst.py
echo "Running: lst -C ex_cal -j 2454555.5 > lst_out.txt"
lst -C ex_cal 2454555.5 | tee lst_out.txt

# Exercise mdlvis.py, create new.uv
echo "Running: mdlvis.py -C ex_cal -s cyg --nchan=16 --sfreq=.1 --sdf=.001 --inttime=60 --startjd=2454555.5 --endjd=2454555.6 --pol=xx -m sim"
mdlvis.py -C ex_cal -s cyg --nchan=16 --sfreq=.1 --sdf=.001 --inttime=60 --startjd=2454555.5 --endjd=2454555.6 --pol=xx -m sim

# Exercise uvlist.py
echo "Running: uvlist.py -k sdf,sfreq,nchan new.uv > uvlist_out.txt"
uvlist.py -k sdf,sfreq,nchan new.uv | tee uvlist_out.txt

# Exercise plot_uv.py
plot_uv.py -p xx -a 0_1 -c 2_14 -x 2 -o plot_uv_out.png new.uv

# Exercise phs2src.py, create new.uv.cyg
echo "Running: phs2src.py -C ex_cal -s cyg"
phs2src.py -C ex_cal -s cyg new.uv

# Exercise combine_freqs.py, create new.uvm
echo "Running: combine_freqs.py -n 8 new.uv"
combine_freqs.py -n 8 new.uv

# Exercise compress_uv.py, create new.uv.tar.gz
echo "Running: compress_uv.py new.uv"
compress_uv.py new.uv

# Exercise xrfi.py, create new.uvr
echo "Running: xrfi.py -c 0_1,14_15 -m none new.uv"
xrfi.py -c 0_1,14_15 -m none new.uv

# Exercise filter_src.py, create new.uvx.cyg
echo "Running: filter_src.py -C ex_cal -s cyg -f 2 -d 2 --clean=1e-3 new.uvm"
filter_src.py -C ex_cal -s cyg -r 2 -d 2 --clean=1e-3 new.uvr

# Exercise difuv.py, create new.uvd
echo "Running: difuv.py new.uv new.uv"
uv_addsub.py --sub new.uv new.uv

# Exercise flux_cal.py, create new.uvf
echo "Running: flux_cal.py -C ex_cal -s cyg -p -b -f new.uv"
flux_cal.py -C ex_cal -s cyg -p -b -f new.uv

# Exercise fitmdl.py, create fitmdl_out.txt
echo 'Running: fitmdl.py -a 0_1 -p xx -c 5_10 -C ex_cal -s cyg -x 10 -P "(0/1)=phsoff" -S "(0/1)=amp" --maxiter=200 new.uv'
fitmdl.py -a 0_1 -p xx -c 5_10 -C ex_cal -s cyg -x 10 -P "(0/1)=phsoff" -S "(0/1)=amp" --maxiter=200 new.uv

# Exercise mk_img.py, create out0.{dim,dbm}.fits
echo "Running: mk_img.py -p xx -C ex_cal -s cyg -o dim,dbm --fmt=out --altmin=30 new.uv"
mk_img.py -p xx -C ex_cal -s cyg -o dim,dbm --fmt=out%d --altmin=30 new.uv

# Exercise cl_img.py, create out0.bim.fits
echo "Running: cl_img.py -d cln -o bim --maxiter=1000 out0.{dim,dbm}.fits"
cl_img.py -d cln -o bim --maxiter=1000 out0.{dim,dbm}.fits

# Exercise plot_img.py, create plot_img_out.png
echo "Running: plot_img.py -o plot_img_out.png out0.bim.fits"
plot_img.py -o plot_img_out.png out0.bim.fits

# Exercise mk_map.py, create mk_map_out.fits
echo "Running: mk_map.py -m mk_map_out.fits --nside=128 out0.bim.fits"
mk_map.py -m mk_map_out.fits --nside=128 out0.bim.fits

# Exercise modmap.py, create modmap_out.fits
echo "Running: modmap.py -C ex_cal -s cas,crab,vir -m modmap_out.fits -i mk_map_out.fits --dtype=float"
modmap.py -C ex_cal -s cas,crab,vir -m modmap_out.fits -i mk_map_out.fits --dtype=float

# Exercise plot_map.py, create plot_map_out.png
echo "Running: plot_map.py -s cyg,cas,crab,vir --drng=4 -p moll -o plot_map_out.png modmap_out.fits"
plot_map.py -s cyg,cas,crab,vir --drng=4 -p moll -o plot_map_out.png modmap_out.fits

# Clean up
echo "Cleaning up..."
rm -rf new.uv* *_out.{txt,png,fits} out0.*.fits

echo "Done."
