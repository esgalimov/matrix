obj="./matrix"

test_folder="../tests/det_test/"

correct=("1" "34" "306" "-9" "24" "2592" "1970.64" "128" "-864" "1024" "7776" "-11664" "-3456")

echo "TESTS:"
echo
for ((i = 1; i <= 13; i++)) do
    echo $i.dat
    ans=`${obj} < ${test_folder}$i.dat`
    echo $ans

    if [ $ans = ${correct[$i-1]} ]
    then
    echo -e "\033[0;32mCORRECT\033[0m"
    else
    echo -e "\033[0;31mERROR\033[0m correct = ${correct[$i-1]}"
    fi
    echo
done
