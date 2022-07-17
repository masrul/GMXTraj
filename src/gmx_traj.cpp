#include "gmx_traj.hpp"
#include <algorithm>

static bool check_file_extension(std::string file_name,std::string ext) {
    size_t i = file_name.rfind('.', file_name.length());
    if (i != std::string::npos) {
        std::string ext_here=file_name.substr(i+1, file_name.length() - i);
        if (ext_here==ext)
            return true;
        else
            return false;
    }
    else 
        return false;
}

static std::string trim(std::string s) {
    // Remove white space from string 
    
    std::string s1=s;   
    s1.erase(s1.begin(), std::find_if_not(s1.begin(), s1.end(), [](char c){ return std::isspace(c); }));
    s1.erase(std::find_if_not(s1.rbegin(), s1.rend(), [](char c){ return std::isspace(c); }).base(), s1.end());
    return s1;
}


GMXTraj::GMXTraj(std::string traj_file,std::string top_file=""){

    // gro as trajectory file and  top_file  is not required
    if (top_file == ""){
        if (check_file_extension(traj_file,"gro")){
            traj_file_extension = "gro";
            _init_top(traj_file);
            _init_gro(traj_file);
        }
        else {
            std::cout << "Trajectory extension must be .gro\n";
            std::terminate();
        }
    }

    // XTC/TRR as traj_file, and top_file is required 
    else {

        if (check_file_extension(top_file,"gro"))
            _init_top(top_file);
        else {
            std::cout << "Topology extension must be .gro\n";
            std::terminate();
        }

        if (check_file_extension(traj_file,"xtc")){
            traj_file_extension = "xtc";
            _init_xtc(traj_file);
        }
        else if (check_file_extension(traj_file,"trr")){
            traj_file_extension = "trr";
            _init_trr(traj_file);
        }
        else {
            std::cout << "Trajectory  extension must be .xtc or .trr\n";
            std::terminate();
        }
    }
}


void GMXTraj::_init_xtc(std::string xtc_file){
    // Check natoms in XTC and initialze XTC
    prec = 1000;
    int natoms_xtc;
    xd = xdrfile_open(&xtc_file[0],"r");
    xdr_status= read_xtc_natoms(&xtc_file[0], &natoms_xtc);

    if (xdr_status == exdrOK){
        if (natoms_xtc !=natoms){
            std::cout << "Number of atoms mismatch between TOP("
                <<natoms<<") and XTC("<<natoms_xtc<<")!\n";
            std::terminate();
        }
        std::cout << "natoms: "<<natoms<<"\n";
        pos = (rvec*)malloc(natoms*3*sizeof(float));
    }
}


void GMXTraj::_init_top(std::string top_file){
    // Colllect symbols, resids,resnames,natoms from GRO file 
    gro_handler.open(top_file);
    if (!gro_handler) {
        std::cout << top_file << " doesn't exist\n";
        std::terminate();
    }
    std::string line;

    getline(gro_handler, line);
    getline(gro_handler, line);
    natoms = stoi(line);
    for (int i = 0; i < this->natoms; ++i) {
        getline(gro_handler, line);
        resnames.push_back(trim(line.substr(5,5)));
        symbols.push_back(trim(line.substr(10,5)));
        resids.push_back(std::stoi(line.substr(0,5)));
    }
    gro_handler.close();
}

void GMXTraj::_init_gro(std::string gro_file){
    gro_handler.open(gro_file);
    pos = (rvec*)malloc(natoms*3*sizeof(float));

    // Set zeros 
    for (int i=0;i<3;++i)
        for (int j=0;j<3;++j)
            box[i][j] = 0.0;
}


void GMXTraj::_init_trr(std::string trr_file){
    // Check natoms in TRR and initialze TRR
    prec = 1000;
    int natoms_trr;
    xd = xdrfile_open(&trr_file[0],"r");
    xdr_status= read_trr_natoms(&trr_file[0], &natoms_trr);

    if (xdr_status == exdrOK){
        if (natoms_trr !=natoms){
            std::cout << "Number of atoms mismatch between TOP("
                <<natoms<<") and TRR("<<natoms_trr<<")!\n";
            std::terminate();
        }
        std::cout << "natoms: "<<natoms<<"\n";
        pos = (rvec*)malloc(natoms*3*sizeof(float));
    }
}


bool GMXTraj::next(){
    int read_status;

    if (traj_file_extension == "xtc")
        read_status=_read_next_xtc();
    else if (traj_file_extension == "trr")
        read_status= _read_next_trr();
    else 
        read_status= _read_next_gro();

    return read_status==READ_SUCCESS ? true : false;
}

int GMXTraj::_read_next_xtc(){
    xdr_status= read_xtc(xd, natoms, &step, &time, box, pos, &prec);
    if (xdr_status==exdrOK)
        return READ_SUCCESS;
    else
        return READ_FAILED;
}

int GMXTraj::_read_next_trr(){
    float lambda =0.0; // ignore 
    int has_prop = 0; // ignore  

    xdr_status = read_trr(
            xd, natoms, &step, &time,&lambda, 
            box, pos,nullptr,nullptr,&has_prop
            );
    if (xdr_status==exdrOK)
        return READ_SUCCESS;
    else
        return READ_FAILED;
    return 0;
}

int GMXTraj::_read_next_gro(){
    if (gro_handler.peek() == EOF)
        return READ_FAILED;

    std::string line;

    // get time
    getline(gro_handler, line);
    size_t pos1 = line.find("t=");
    size_t pos2 = line.find("step=");
    std::string target1 = line.substr(pos1 + 2, pos2 - pos1 - 2);
    if (pos1 ==std::string::npos || pos2 ==std::string::npos)
        time = 0.0;
    else 
        time = stof(target1);
    step = 0 ; // need to implement 

    getline(gro_handler, line);

    for (int i = 0; i < natoms; ++i) {
        getline(gro_handler, line);
        pos[i][0] = stof(line.substr(20,8));
        pos[i][1] = stof(line.substr(28,8));
        pos[i][2] = stof(line.substr(36,8));
    }

    getline(gro_handler, line);
    std::istringstream input(line);
    input >> lx >> ly >> lz;

    return READ_SUCCESS;
}

void GMXTraj::create_residue_tracker(){
    residue_trackers.clear();

    for (int i=0; i< natoms; ++i){
        int resid = resids[i];
        std::string resname = resnames[i];

        int rIDx = residue_trackers.size()-1;

        if (residue_trackers.size() == 0){
            Tracker residue_tracker(resid,resname,1,i);
            residue_trackers.push_back(residue_tracker);
        }
        else if 
            ((residue_trackers[rIDx].id == resid) & 
            (residue_trackers[rIDx].name == resname)){
            
            residue_trackers[rIDx].natoms +=1;
        }
        else {
            Tracker residue_tracker(resid,resname,1,i);
            residue_trackers.push_back(residue_tracker);
        }
    }
}




void GMXTraj::create_molecule_tracker(std::string sysinfo_file=""){
    if (sysinfo_file == "")
        _molecule_tracker_from_residue();
    else 
        _molecule_tracker_from_file(sysinfo_file);
}


void GMXTraj::_molecule_tracker_from_residue(){
    // Works if each molecule has exactly one residue 
    molecule_trackers.clear();

    for (int i=0; i< natoms; ++i){
        int molid = resids[i];
        std::string molname = resnames[i];

        int mIDx = molecule_trackers.size()-1;

        if (molecule_trackers.size() == 0){
            Tracker molecule_tracker(molid,molname,1,i);
            molecule_trackers.push_back(molecule_tracker);
        }
        else if 
            ((molecule_trackers[mIDx].id == molid) & 
            (molecule_trackers[mIDx].name == molname)){
            
            molecule_trackers[mIDx].natoms +=1;
        }
        else {
            Tracker molecule_tracker(molid,molname,1,i);
            molecule_trackers.push_back(molecule_tracker);
        }
    }
    nmolecules = molecule_trackers.size();
}

void GMXTraj::_molecule_tracker_from_file(std::string sysinfo_file){
    molecule_trackers.clear();
    std::ifstream file; 
    file.open(sysinfo_file);
    if (!file) {
        std::cout << sysinfo_file << " doesn't exist\n";
        std::terminate();
    }

    std::string line;
    std::vector<MoleculeSummary> molecule_summarys;

    while(getline(file,line)){
        if (line !=""){
            std::istringstream input(line);
            MoleculeSummary summary; 
            input >> summary.name >> summary.natoms >> summary.nitems;
            molecule_summarys.push_back(summary);
        }
    }
    file.close();

    int id =1; 
    int sIDx = 0;
    int _natoms = 0;
    
    for (auto summary : molecule_summarys){
        for (int i=0; i< summary.nitems; ++i){
            molecule_trackers.push_back(Tracker(id,summary.name,summary.natoms,sIDx));
            sIDx +=summary.natoms;
            id +=1;
            _natoms +=summary.natoms;
        }
    }
    assert(_natoms==natoms && "Number of atoms mismatch between topology and molecule summary");
    nmolecules = molecule_trackers.size();
}


GMXTraj::~GMXTraj(){

    if (pos != nullptr)
        delete [] pos;
}


