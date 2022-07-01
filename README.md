# GMXTraj

**GMXTraj** is C++ library that wraps *xdrfile* with object oriented fahsion. 

```cpp
#include "gmx_traj.hpp"
#include <string>

int main(){
    std::string top_file; // gro file 
    std::string traj_file; // xtc, trr file 
    traj = GMXTraj(top_file,traj_file); 

    // Iterate trajectory
    while (traj.next()){
        // do staff 
        // traj.pos to access position 
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
```
cd GMXTraj 
mkdir build  && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install 
make 
make install 
```

+ Install path contains header files, add **gmx_traj.hpp** into your *cpp* file 

+ Best way to link library using  **CMakeLists.txt** contains in example folder. Edit 
according to your needs 

