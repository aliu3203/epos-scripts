cd /media/Disk_Jin/aliu
cd epos1
i="$1"
mv "z-auau_run_1.root" "z-auau_run_${i}.root"
mv "z-auau_run_${i}.root" ../auau200-7
cd ../epos2
i=$((i+1))
mv "z-auau_run_2.root" "z-auau_run_${i}.root"
mv "z-auau_run_${i}.root" ../auau200-7
cd ../epos3
i=$((i+1))
mv "z-auau_run_3.root" "z-auau_run_${i}.root"
mv "z-auau_run_${i}.root" ../auau200-7
cd ../epos4
i=$((i+1))
mv "z-auau_run_4.root" "z-auau_run_${i}.root"
mv "z-auau_run_${i}.root" ../auau200-7
cd ../epos5
i=$((i+1))
mv "z-auau_run_5.root" "z-auau_run_${i}.root"
mv "z-auau_run_${i}.root" ../auau200-7
