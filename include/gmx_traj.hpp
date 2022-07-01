#ifndef GMX_TRAJ_H

#define GMX_TRAJ_H

#include<iostream>
#include<vector>
#include <string>
#include<fstream>
#include <sstream> 
#include<iostream>
#include<vector>
#include <string>
#include<cassert>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"
#include <limits>

struct Tracker{
    int id;
    std::string name;
    int natoms; // natoms per residue or molecule 
    int sIDx; // strating index of this residue or molecule 
    Tracker(int,std::string,int,int);
};
Tracker::Tracker(int _id,std::string _name,int _natoms,int _sIDx){
    id = _id;
    name = _name;
    natoms= _natoms;
    sIDx = _sIDx;
}
typedef struct Tracker Tracker;



struct MoleculeSummary{
    std::string name;
    int natoms; 
    int nitems;
    MoleculeSummary(std::string,int,int);
};
MoleculeSummary::MoleculeSummary(std::string _name,int _natoms,int _nitems){
    name = _name;
    natoms = _natoms;
    nitems = _nitems;
}
typedef struct MoleculeSummary MoleculeSummary; 



enum {READ_SUCCESS,READ_FAILED};

class GMXTraj{
public:
    int natoms;
    int nmolecules;
    rvec* pos;
    float time;
    int step;
    matrix box;
    float &lx=box[0][0]; 
    float &ly=box[1][1]; 
    float &lz=box[2][2]; 
    std::vector<std::string> symbols;
    std::vector<std::string> resnames;
    std::vector<int> resids;
    std::vector<Tracker> residue_trackers;
    std::vector<Tracker> molecule_trackers;
    
    // Constructors 
    GMXTraj(std::string,std::string);
    GMXTraj(std::string);

    // Destructors
    ~GMXTraj();

    // Followings should be called by users 
    bool next();  // get next frame 
    void create_residue_tracker();
    void create_molecule_tracker(std::vector<MoleculeSummary>);
    void make_whole(std::string); // fix borken residue or molecules; 

private:
    int xdr_status;
    float prec;
    std::string traj_file_extension; 
    std::ifstream gro_handler;
    XDRFILE* xd;
    void _init_top(std::string);
    void _init_xtc(std::string);
    void _init_trr(std::string);
    void _init_gro(std::string);
    int _read_next_xtc();
    int _read_next_trr();
    int _read_next_gro();
    void _make_whole_residue();
    void _make_whole_molecule();
};



#endif /* end of include guard GMX_TRAJ_H */

