# GMXTraj

**GMXTraj** is a C++ library that wraps *xdrfile* with object oriented fahsion. GMXTraj can be used 
mainly for analysis purpose. 

```cpp
#include "gmx_traj.hpp"
#include <string>

int main(){
    std::string top_file; // gro file 
    std::string traj_file; // gro,xtc, trr file 
    GMXTraj traj{traj_file,top_file}; // traj{traj_file} also vaild for gro  

    // Iterate trajectory
    while (traj.next()){
        // do staff 
        // traj.pos to access position e.g. traj.pos[5][0] for x of 6th atom 
        // traj.box to access box matrix
        // traj.lx, traj.ly, traj.lz can be used for orthagonal box lengths
    }
    return 0;
}
```

## How to use 

+ Clone directory 
```bash 
git clone https://github.com/masrul/GMXTraj
```
+ Build library 
```bash 
cd GMXTraj 
mkdir build  && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install 
make 
make install 
```

+ Install path contains header files, add **gmx_traj.hpp** into your *cpp* file 

+ Best way to link library using  **CMakeLists.txt** located at example folder. Edit 
according to your needs 

