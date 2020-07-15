#include "correlation.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "VectorMath.h"

std::vector<double> orient_corr(const std::vector<std::vector<double> >& angles);

std::vector<double> pos_orient_corr(const std::vector<std::vector<vec2<double> > >&, const std::vector<std::vector<double> >&, const double&);

int main(int argc, char** argv){
    std::string filename = (argc > 1)? argv[1] :"./Result.txt";
    std::ifstream input( filename);
    bool load_vector = false;
    int numberParticles;
    double boxSize;
    std::string line;
    std::getline(input, line);
    std::istringstream sline(line);
    std::string tempstring;
    std::getline(sline, tempstring, ';');
    numberParticles = stoi(tempstring);
    std::getline(sline, tempstring, ';');
    boxSize = stod(tempstring);
    //std::cout<< boxSize << std::endl;


    std::vector<std::vector<vec2<double> > > positions(numberParticles, std::vector<vec2<double> >());
    std::vector<std::vector<double> > orientations(numberParticles, std::vector<double>());
    int idx = 0;
    while(input){
        std::string line;
        getline(input, line);
        std::istringstream sline(line);
        int col = 0;
        int curr_idx = idx % numberParticles;
        double previous;
        while(sline){
            std::string tempstring;
            if(!std::getline(sline, tempstring, ';')) break;
            double number = std::stod(tempstring);
            if(col == 1){
                positions[curr_idx].push_back(vec2<double>(previous, number));
            }
            else if(col == 2){
                orientations[curr_idx].push_back(number);
            }
            previous = number;
            ++col;
        }
        ++idx;
    }
    input.close();

    /*for(int i = 0; i < 100; ++i){
        std::cout << positions[0][i].x << std::endl;
        std::cout << orientations[0][i] << std::endl;
    }*/
    //std::cout<<dot(vec2<double>(positions[0][1] - positions[0][0]) ,vec2<double>(cos(orientations[0][0]), sin(orientations[0][0]))) << std::endl;
    std::string orOutput = filename.substr(0, filename.size()-4) + "_orCorr.txt";
    std::ofstream output(orOutput);
    std::vector<double> orCorr = orient_corr(orientations);
    for(auto cor: orCorr){
        output<< cor << "\n";
    }
    output.close();
    std::string porOutput = filename.substr(0, filename.size()-4) + "_porOrCorr.txt";
    std::ofstream outfile(porOutput);
    std::vector<double> posOrCorr = pos_orient_corr(positions, orientations, boxSize);
    for(auto cor:posOrCorr){
        outfile << cor << "\n";
    }
    outfile.close();

}


std::vector<double> orient_corr(const std::vector<std::vector<double> >& angles){
    const unsigned int numberParticles = angles.size();
    const unsigned int timeLength = angles[0].size();
    std::vector<double> res(timeLength/2, 0);
    for(int idx = 0; idx < numberParticles; ++idx){
        for(int deltaT = 0; deltaT < timeLength/2; ++deltaT){
            double cor_sum = 0, var = 0;
            unsigned int range = timeLength - deltaT;
            for(int t = 0; t < range; ++t){
                cor_sum += cos(angles[idx][t+deltaT] - angles[idx][t]);
                var += 1;
            }
            res[deltaT] += cor_sum/var;  
        }
    }
    for(int deltaT = 0; deltaT < timeLength/2; ++deltaT){
        res[deltaT] /= numberParticles;
    }
    return res;
}


std::vector<double> pos_orient_corr(const std::vector<std::vector<vec2<double> > >& positions, const std::vector<std::vector<double> >& orientations, const double& boxSize){
    int numParticles = orientations.size();
    int timeLength = orientations[0].size();
    int deltaLength = timeLength/2;
    std::vector<double> res(deltaLength, 0);
    int cross_boundary = 0;
    for(int idx = 0; idx < numParticles; ++idx){
        /*
        for(int deltaT = 0; deltaT < deltaLength; ++ deltaT){
            double cor_sum = 0;
            for(int t = 0; t < deltaLength; ++ t){
                vec2<double> disp(positions[idx][t+deltaT] - positions[idx][t]);
                vec2<double> angle(cos(orientations[idx][t]), sin(orientations[idx][t]));
                cor_sum += dot(disp, angle);
            }
            res[deltaT] += cor_sum/deltaLength;
        }*/
        for(int t = 0; t < deltaLength; ++t){
            //double cor_sum = 0, var = 0;
            //unsigned int range = timeLength - deltaT;
            vec2<double> previous;
            for(int deltaT = 0; deltaT < deltaLength; ++deltaT){    
                vec2<double> disp(positions[idx][t+deltaT] - positions[idx][t]);
                if(dot(disp, previous) < 0){
                    vec2<double> velocity(disp - previous);
                    if(velocity.x < -0.5*boxSize) {disp.x += boxSize;}
                    else if(velocity.x > 0.5*boxSize) disp.x -= boxSize;
                    if(velocity.y < -0.5*boxSize) disp.y += boxSize;
                    else if(velocity.y > 0.5*boxSize) disp.y -= boxSize;
                }
                previous = disp;
                vec2<double> angle(cos(orientations[idx][t]), sin(orientations[idx][t]));
                res[deltaT] += dot(disp, angle);
            }
        }
    }
    int norm = (numParticles*deltaLength);
    for(int deltaT = 0; deltaT < res.size(); ++deltaT){
        res[deltaT] /= norm;
    }
    //std::cout << cross_boundary << std::endl;
    return res;
}

