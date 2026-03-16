set debuginfod enabled off
set confirm off

file ./build/tests/edt
set args example.db

break distance.cpp:61
condition 1 i== 5 * 10*10 + 8 * 10 +  5

run
