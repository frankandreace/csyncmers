echo "Disabling HYPERTREADING for more accurate time estimation"
sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'

klen=31
slen=11
filename="data/chr19_bit.fa"
# "data/chr19_bit.fa"
echo "./bin/syncmer \"$filename\" $klen $slen"
./bin/syncmer $filename $klen $slen

echo "Re-enabling HYPERTREADING."
sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'