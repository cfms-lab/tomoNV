# TomoNV_Win64_DLL
### C++ implementation of the tomoNV algorithm


- tested in Visual Studio Community 2019
- To run the DLL, need "Cpp_example.py" of https://github.com/cfms-lab/tomoNV.
- Compile and copy the generated "TomoNV_Win64.dll" and "TomoNV_Win64.lib" files to the folder where "Cpp_example.py" is.
- Number of mesh triangles should be smaller than 2,147,483,647 (Using int32 for triangle indices in  "Tomo_types.h")
```
DLLEXPORT typedef int         MESH_ELE_ID_TYPE;
```
- The interface function to Python is "TomoNV_TMPxl()" of "TomoNV_Win64.h".
- To make a standalone application(.exe) based on these souce codes, delete the definition "\_CREATING_DLL\_" in "Configuration Properties->C/C++->Preprocessor->Preprocessor Definitions" of Project Properties window.
