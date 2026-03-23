set debuginfod enabled off
set confirm off

file ./build/tests/edt
set args sphere-pack-50c-uint8.db

set $X = 1
set $Y = 49
set $Z = 48

break distance.cpp:61
condition 1 i == $X + Distance.size(0) * ($Y+ $Z * Distance.size(1) )
run
