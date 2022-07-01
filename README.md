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
