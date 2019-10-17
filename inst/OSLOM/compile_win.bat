

set Main_folder=Sources_2_5

set source_folder=%Main_folder%\OSLOM_files

g++ -static -o oslom_undir_win %source_folder%\main_undirected.cpp -O3 -Wall


