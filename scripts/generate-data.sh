# Launch with arguments "model orders"
model=${1-radial}
order=${2-"3,5,7,9"}

set -f          # disable globbing

# extract all results into two files (distortion/correction)
rm -f $model.dat ${model}r.dat
IFS=$'\n' # make newlines the only separator
for i in $(cat models.txt);
do
   echo \"$i\"
   ../precision_analysis -d \"$i\" $model $order >> $model.dat
   ../precision_analysis -d \"$i\" $model $order -r >> ${model}r.dat
done
IFS=$' \t\n'

# demultiplex RMSE (field 3) into different files, one per order, sort them
for i in $(echo $order | tr , \ );
do
  grep "^$i " $model.dat    | cut -d \  -f 3 > $model$i.dat
  grep "^$i " ${model}r.dat | cut -d \  -f 3 > $model${i}r.dat
  sort -g $model$i.dat    > $model${i}s.dat
  sort -g $model${i}r.dat > $model${i}rs.dat
done

# demultiplex max (field 4) into different files, one per order, sort them
for i in $(echo $order | tr , \ );
do
  grep "^$i " $model.dat    | cut -d \  -f 4 > $model${i}m.dat
  grep "^$i " ${model}r.dat | cut -d \  -f 4 > $model${i}rm.dat
  sort -g $model${i}m.dat  > $model${i}ms.dat
  sort -g $model${i}rm.dat > $model${i}rms.dat
done
