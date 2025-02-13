cmake_minimum_required(VERSION 3.5)

project(quartet LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(quartet 
 	src/predicates.cpp           
	src/geometry_queries.cpp     
	src/sdf.cpp                  
	src/trimesh.cpp              
	src/tet_mesh.cpp             
	src/feature.cpp              
	src/read_obj.cpp             
	src/tet_quality.cpp          
	src/match_features.cpp       
	src/optimize_tet_mesh.cpp    
	src/make_signed_distance.cpp 
	src/make_tet_mesh.cpp        
	src/main.cpp
)

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
	target_link_libraries(quartet PUBLIC OpenMP::OpenMP_CXX)
endif()

# fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp omptest.c -o omptest 

