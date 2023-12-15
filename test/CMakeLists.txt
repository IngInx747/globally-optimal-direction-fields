# main
add_executable("ffgen" "main.cc")
target_link_libraries("ffgen" libffgen)
add_dependencies("ffgen" libffgen)
