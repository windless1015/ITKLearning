cmake_minimum_required(VERSION 3.15)
project(itkImagePointer)
set(ITK_IMAGE_POINTER ${PROJECT_NAME})

set(VTK_DIR D:/projects/ultras/libs/VTK/x64/lib/cmake/vtk-8.2)
set(ITK_DIR D:/projects/ultras/libs/ITK/x64/lib/cmake/ITK-5.1)

find_package(ITK)
find_package(VTK REQUIRED)


#set(ITK_REQUIRED_MODULES "ITKCommon;ITKIOImageBase;ITKIOJPEG;ITKIOGDCM;ITKCommon;ITKImageIntensity;ITKVTK;ITKVtkGlue;ITKIONRRD")

include_directories(
${ITK_INCLUDE_DIRS}
${VTK_INCLUDE_DIRS}
)


add_executable(${ITK_IMAGE_POINTER} itkImagePointer.cpp)
target_link_libraries(${ITK_IMAGE_POINTER} 
 ${ITK_LIBRARIES}
 ${VTK_LIBRARIES}
)