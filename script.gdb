set debuginfod enabled off
set confirm off

file ./build/tests/edt
set args example.db

break distance.cpp:61
condition 1 i == 3 + Distance.size(0) * (1+ 3 * Distance.size(1) )
run
