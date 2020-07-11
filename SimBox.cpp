#include "VectorMath.h"
#include "geometry.h"
#include "aabb/AABB.h"
#include "myrandom.h"
#include <random>
#include <fstream>

using namespace std;


inline void periodBoundary(vector<double>& , vector<bool>& , vector<double>& );

int main(int argc, char** argv){
    const unsigned int nIteration = (argc > 3)? stoi(argv[3]) : 200;
    const unsigned int numParticles = (argc > 1) ? stoi(argv[1]) : 100;
    const unsigned int activeParticles = numParticles;
    const double density = (argc > 2)? stod(argv[2]) : 0.25;
    const unsigned int stepSave = (argc > 4)? stoi(argv[4]) : 1;
    const double maxDisp = 0.1;
    const double sq2 = sqrt(0.5);
    const double PI2 = M_PI * 2;
    const double diagRadius = sqrt(0.5);
    double D0 = 0.092122;
    double DR = 0.238;
    double DT = 0.1;
    double v0 = 0.063;
    double scaleSize = 0.99;
    double sidelength = sqrt(numParticles / density);
    vector<double> boxSize({sidelength, sidelength});
    vec2<double> box(sidelength, sidelength);
    vector<bool> period({true, true});
    unsigned int success_move = 0;
    unsigned int activeCount = 0;
    Rand myGenerator;

    ofstream file("./Result.txt");
    file << numParticles <<";" << boxSize[0] << ";" << boxSize[1] << '\n';
    //output.close();
    /***********************************************************/
    /*                  Initialize System                      */
    /***********************************************************/
    aabb::Tree simTree(2, maxDisp, period, boxSize, numParticles);
    vector<vector<double> > positions;
    vector<double> orientations;
    vector<Square> squareList;
    vector<unsigned int> tags;
    vector<unsigned int> idxs;
    
    vector<double> activeNoises(numParticles, D0);
    vector<double> transverseNoises(numParticles, DT);
    vector<double> rotationNoises(numParticles, DR);
    vector<double> activeVelocity(numParticles, v0);
    vector<int> failedMoves(numParticles, 0);

    for(unsigned int i = 0; i < numParticles; ++i){
        vector<double> pos(2);
        double orientation;
        Square tempSquare;
        if(i == 0){
            pos[0] = myGenerator.rdnm()*boxSize[0];
            pos[1] = myGenerator.rdnm()*boxSize[1];

            orientation = myGenerator.rdnm() * PI2;

            tempSquare = Square(pos, orientation);
        }
        else{
            bool overlap = true;

            while(overlap){
                pos[0] = myGenerator.rdnm()*boxSize[0];
                pos[1] = myGenerator.rdnm()*boxSize[1];

                orientation = myGenerator.rdnm()*PI2;

                tempSquare = Square(pos, orientation);

                vector<double> lowerBound({pos[0] - diagRadius, pos[1] - diagRadius});
                vector<double> upperBound({pos[0] + diagRadius, pos[1] + diagRadius});

                aabb::AABB singleAABB(lowerBound, upperBound);

                vector<unsigned int> particles = simTree.query(singleAABB);

                overlap = false;
                for(unsigned int j = 0; j < particles.size(); ++j){
                    if(test_polygon_overlap(squareList[particles[j]], tempSquare, box)){
                        overlap = true;
                        break;
                    }
                }
            }
        }
        simTree.insertParticle(i, pos, diagRadius);
        orientations.push_back(orientation);
        positions.push_back(pos);
        squareList.push_back(tempSquare);
        unsigned int tag = activeCount < activeParticles ? 1 : 0;
        tags.push_back(tag);
        ++activeCount; 
        idxs.push_back(i);
    }
    //write position and file to file;
    //fstream file;
    //file.open("./Result.txt", fstream::app);
    for(int i = 0; i < positions.size(); ++i){
        file << positions[i][0] << ";" << positions[i][1] << ";" << orientations[i] << "\n";
    }

    /***********************************************************/
    /*                  Starting Dynamics                      */
    /***********************************************************/
    
    for(int it = 0; it < nIteration; ++it){
        random_shuffle(idxs.begin(), idxs.end());
        for(unsigned int idx: idxs){
            unsigned int pt_tag = tags[idx];
            double curr_orient = orientations[idx];

            double mu = (pt_tag == 1) ? activeVelocity[idx]: 0;

            double rotation_disp = myGenerator.normal_rdnm(0, rotationNoises[idx]);
            double disp_0 = myGenerator.normal_rdnm(mu, activeNoises[idx]);
            double disp_t = myGenerator.normal_rdnm(0, transverseNoises[idx]);

            double disp_x = disp_0 * cos(curr_orient) - disp_t * sin(curr_orient);
            double disp_y = disp_0 * sin(curr_orient) + disp_t * cos(curr_orient);

            vector<double> pos({positions[idx][0] + disp_x, positions[idx][1] + disp_y});

            periodBoundary(pos, period, boxSize);

            Square tempSquare(pos, curr_orient + rotation_disp);
            vector<double> lowerBound({pos[0] - diagRadius, pos[1] - diagRadius});
            vector<double> upperBound({pos[0] + diagRadius, pos[1] + diagRadius});
            
            aabb::AABB boundBox(lowerBound, upperBound);

            vector<unsigned int> particles = simTree.query(boundBox);

            bool overlap = false;

            for(unsigned int pt: particles){
                if((pt != idx)&&test_polygon_overlap(squareList[pt], tempSquare, box)){
                    overlap = true;
                    break;
                }
            }

            if(!overlap){
                ++success_move;
                positions[idx] = pos;
                squareList[idx] = tempSquare;
                orientations[idx] = curr_orient + rotation_disp;
                simTree.updateParticle(idx, lowerBound, upperBound);
                failedMoves[idx] = 0; 
            }
            else{
                disp_x = disp_0 * sq2 * cos(curr_orient - M_PI_4) - disp_t * sq2 * sin(curr_orient + M_PI_4);
                disp_y = disp_0 * sq2 * sin(curr_orient - M_PI_4) + disp_t * sq2 * cos(curr_orient + M_PI_4);

                pos = vector<double>({positions[idx][0] + disp_x, positions[idx][1] + disp_y});

                periodBoundary(pos, period, boxSize);

                tempSquare = Square(pos, curr_orient + rotation_disp);
                lowerBound = vector<double>({pos[0] - diagRadius, pos[1] - diagRadius});
                upperBound = vector<double>({pos[0] + diagRadius, pos[1] + diagRadius});
                boundBox = aabb::AABB(lowerBound, upperBound);
                particles = simTree.query(boundBox);
                overlap = false;
                for(unsigned int pt: particles){
                    if((pt != idx)&&test_polygon_overlap(squareList[pt], tempSquare, box)){
                        overlap = true;
                        break;
                    }
                }
                if(!overlap){
                    ++success_move;
                    positions[idx] = pos;
                    squareList[idx] = tempSquare;
                    orientations[idx] = curr_orient + rotation_disp;
                    simTree.updateParticle(idx, lowerBound, upperBound);
                    failedMoves[idx] = 0;
                }
                else{
                    disp_x = disp_0 * sq2 * cos(curr_orient + M_PI_4) - disp_t * sq2 * sin(curr_orient - M_PI_4);
                    disp_y = disp_0 * sq2 * sin(curr_orient + M_PI_4) + disp_t * sq2 * cos(curr_orient - M_PI_4);

                    pos = vector<double>({positions[idx][0] + disp_x, positions[idx][1] + disp_y});

                    periodBoundary(pos, period, boxSize);

                    tempSquare = Square(pos, curr_orient + rotation_disp);
                    lowerBound = vector<double>({pos[0] - diagRadius, pos[1] - diagRadius});
                    upperBound = vector<double>({pos[0] + diagRadius, pos[1] + diagRadius});
                    boundBox = aabb::AABB(lowerBound, upperBound);
                    particles = simTree.query(boundBox);
                    overlap = false;
                    for(unsigned int pt: particles){
                        if((pt != idx)&&test_polygon_overlap(squareList[pt], tempSquare, box)){
                            overlap = true;
                            break;
                        }
                    }
                    if(!overlap){
                        ++success_move;
                        positions[idx] = pos;
                        squareList[idx] = tempSquare;
                        orientations[idx] = curr_orient + rotation_disp;
                        simTree.updateParticle(idx, lowerBound, upperBound);
                        failedMoves[idx] = 0;
                    }
                    else{
                        ++failedMoves[idx];
                        activeNoises[idx]*=scaleSize;
                        transverseNoises[idx]*=scaleSize;
                        rotationNoises[idx]*=scaleSize;
                        activeVelocity[idx]*=scaleSize;
                    }                  
                }
            if(!failedMoves[idx] && (activeNoises[0] < D0)){
                activeNoises[idx]/=scaleSize;
                transverseNoises[idx]/=scaleSize;
                rotationNoises[idx]/=scaleSize;
                activeVelocity[idx]/=scaleSize;
            }

            }
        }
        std::cout << "Finish Iteration " << (it+1) << endl;
        if(it % stepSave == 0){
            std::cout << it << endl;
            for(int i = 0; i < numParticles; ++i){
                file << positions[i][0] << ";" << positions[i][1] << ";" << orientations[i] << "\n";
            }
        }

    }
    file.close();
}






inline void periodBoundary(vector<double>& pos, vector<bool>& period, vector<double>& boxSize){
    for(int i = 0; i < 2; ++i){
        if(pos[i] < 0) pos[i] += period[i] * boxSize[i];
        else if(pos[i] >= boxSize[i]) pos[i] -= period[i]* boxSize[i];
    }
    return;
}


