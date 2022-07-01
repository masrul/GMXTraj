/**
 * @author      : Masrul Huda (mail2masrul@gmail.com)
 * @file        : main
 * @created     : Thursday Jun 30, 2022 11:22:35 CDT
 */
#include "gmx_traj.hpp"

int main(){

    /* { */
    GMXTraj traj("../nvt.xtc","../nvt.gro");

    while (traj.next()){
        std::cout << "Steps: "<<traj.step <<  " Time: "<<traj.time<<"\n";
        
        for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
                std::cout << traj.box[i][j] << " ";
        std::cout << "\n";
    }
    std::cout << traj.lx << " "<<traj.ly<<" "<<traj.lz<<"\n";
    /* } */

    /* { */
    /* GMXTraj traj("nvt.trr","nvt.gro"); */

    /* while (traj.next()==READ_SUCCESS){ */
    /*     std::cout << "Steps: "<<traj.step <<  " Time: "<<traj.time<<"\n"; */
        
    /*     for (int i=0;i<3;++i) */
    /*         for (int j=0;j<3;++j) */
    /*             std::cout << traj.box[i][j] << " "; */
    /*     std::cout << "\n"; */
    /* } */
    /* std::cout << traj.lx << " "<<traj.ly<<" "<<traj.lz<<"\n"; */
    /* } */

    /* GMXTraj traj("nvt_trj.gro"); */

    /* while (traj.next()==READ_SUCCESS){ */
    /*     std::cout << "Steps: "<<traj.step <<  " Time: "<<traj.time<<"\n"; */
        
    /*     for (int i=0;i<3;++i) */
    /*         for (int j=0;j<3;++j) */
    /*             std::cout << traj.box[i][j] << " "; */
    /*     std::cout << "\n"; */
    /* } */
    /* std::cout << traj.lx << " "<<traj.ly<<" "<<traj.lz<<"\n"; */
    
    /* traj.create_residue_tracker(); */

    /* for (int i=0; i<traj.residue_trackers.size();++i){ */
    /*     std::cout */ 
    /*         << "Name: "<< traj.residue_trackers[i].name << " " */
    /*         << "ID: "<< traj.residue_trackers[i].id << " " */
    /*         << "nAtoms: "<< traj.residue_trackers[i].natoms << " " */
    /*         << "sIDx: "<< traj.residue_trackers[i].sIDx << "\n"; */
    /* } */

    /* std::vector<MoleculeSummary> molecule_summarys; */
    /* molecule_summarys.push_back(MoleculeSummary("MFI",72336,1)); */
    /* molecule_summarys.push_back(MoleculeSummary("Lignin",1513,1)); */
    /* molecule_summarys.push_back(MoleculeSummary("Methanol",6,30500)); */
    /* molecule_summarys.push_back(MoleculeSummary("Water",4,30500)); */
    /* traj.create_molecule_tracker(molecule_summarys); */


    /* for (int i=0; i<traj.molecule_trackers.size();++i){ */
    /*     std::cout */ 
    /*         << "Name: "<< traj.molecule_trackers[i].name << " " */
    /*         << "ID: "<< traj.molecule_trackers[i].id << " " */
    /*         << "nAtoms: "<< traj.molecule_trackers[i].natoms << " " */
    /*         << "sIDx: "<< traj.molecule_trackers[i].sIDx << "\n"; */
    /* } */

    return 0;
}


