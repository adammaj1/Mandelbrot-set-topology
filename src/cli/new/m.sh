#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as d.sh
# chmod +x d.sh
# ./d.sh
# checked in https://www.shellcheck.net/




printf "make pgm files \n"
gcc m.c -lm -Wall -Wextra -march=native -fopenmp

if [ $? -ne 0 ]
then
    echo ERROR: compilation failed !!!!!!
    exit 1
fi


export  OMP_DISPLAY_ENV="TRUE"
printf "display OMP info \n"

# https://stackoverflow.com/questions/7526971/how-to-redirect-both-stdout-and-stderr-to-a-file
printf "run the compiled program\n"
time ./a.out > m.txt # 2>> error.txt   # >> adds

export  OMP_DISPLAY_ENV="FALSE"

printf "change Image Magic settings\n"
export MAGICK_WIDTH_LIMIT=100MP
export MAGICK_HEIGHT_LIMIT=100MP

printf "convert all pgm files to png using Image Magic v 6 convert \n"
# for all pgm files in this directory
for file in *.pgm ; do
  # b is name of file without extension
  b=$(basename "$file" .pgm)
  # convert  using ImageMagic : -resize widthxheight || 
  convert "${b}".pgm -resize 600x400 "${b}".png  # iWidth = iHeight* DisplayAspectRatio -resize '10%'
  echo "$file"
done


printf "delete all pgm files \n"
rm ./*.pgm

 
echo OK

printf "info about software \n"
bash --version
make -v
gcc --version
convert -version
convert -list resource
# end
