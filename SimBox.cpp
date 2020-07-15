#include "VectorMath.h"
#include "geometry.h"
#include "aabb/AABB.h"
#include "myrandom.h"
#include <random>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cxxopts.hpp>

using namespace std;


inline void periodBoundary(vector<double>& , vector<bool>& , vector<double>& );

int main(int argc, char** argv){

    cxxopts::Options options("Square Sim", "Simulation of Active Squares in 2D System");

    options.add_options()
        ("n, nIteration", "number of iteration", cxxopts::value<int>()->default_value("10000"))
        ("p, nParticle", "number of total particles", cxxopts::value<int>()->default_value("1000"))
        ("r, activeRatio", "ratio of active particles", cxxopts::value<double>()->default_value("1"))
        ("d, density", "density of system", cxxopts::value<double>()->default_value("0.4"))
        ("t, totalConfiguration", "total configuration to save", cxxopts::value<int>()->default_value("200"))
        ("d0", "diffusion along active direction", cxxopts::value<double>()->default_value("0.00424"))
        ("dt", "transverse diffusion constant", cxxopts::value<double>()->default_value("0.005"))
        ("dr", "rotation diffusion noise", cxxopts::value<double>()->default_value("0.0283"))
        ("v0", "self propelled velocity", cxxopts::value<double>()->default_value("0.063"))
        ("c, continue", "continue previous simulation",  cxxopts::value<bool>()->default_value("false"))
        ("l, loadPath", "load previous result path", cxxopts::value<string>())
    ;

    auto parseResult = options.parse(argc, argv);

    //const unsigned int nIteration = (argc > 3)? stoi(argv[3]) : 200;
    const unsigned int nIteration = parseResult["nIteration"].as<int>();
    unsigned int numParticles = parseResult["nParticle"].as<int>();
    double ratio = parseResult["activeRatio"].as<double>();
    const unsigned int activeParticles = numParticles*ratio;
    const double density = parseResult["density"].as<double>();
    const unsigned int stepSave = nIteration/parseResult["totalConfiguration"].as<int>();
    const double maxDisp = 0.1;
    const double sq2 = sqrt(0.5);
    const double PI2 = M_PI * 2;
    const double diagRadius = sqrt(0.5);
    double real_D0 =  parseResult["d0"].as<double>();
    double D0 = sqrt(2 * real_D0);
    double real_DR = parseResult["dr"].as<double>();
    double DR = sqrt(2 * real_DR);
    double real_DT = parseResult["dt"].as<double>();
    double DT = sqrt(2 * real_DT);
    double v0 = parseResult["v0"].as<double>();
    double sidelength = sqrt(numParticles / density);
    vector<double> boxSize({sidelength, sidelength});
    vec2<double> box(sidelength, sidelength);
    vector<bool> period({true, true});
    unsigned int success_move = 0;
    unsigned int activeCount = 0;
    Rand myGenerator;
    bool loadPrevious = parseResult["continue"].as<bool>();
    string loadPath = (loadPrevious) ? parseResult["loadPath"].as<string>() : "";

    /***********************************************************/
    /*                 Change Name of File                     */
    /***********************************************************/    
    
    stringstream stream;

    double perLen = 2*v0/(DR*DR);
    cout << perLen << endl;

    string filename = "./Result/Result_";

    vector<double> names({perLen, density, v0, real_DR, real_DT});
    for(vector<double>::iterator it = names.begin(); it != names.end(); ++it){
        stream.str(string());
        stream << std::fixed << std::setprecision(3) << *it;
        string s = stream.str();
        filename += (s + "_");
    }
    filename += (to_string(nIteration) + ".txt");


    /***********************************************************/
    /*                     Save Base Info                      */
    /***********************************************************/   

    fstream file;
    if(!loadPrevious){
        file.open(filename, fstream::out);
        file << numParticles <<";" << boxSize[0] << ";" << boxSize[1] << '\n';
    }
    else{
        cout << loadPath << endl;
        file.open(loadPath, fstream::in | fstream::app);
        string line;
        getline(file, line);
        std::istringstream sline(line);
        std::string tempstring;
        std::getline(sline, tempstring, ';');
        numParticles = stoi(tempstring);
        std::getline(sline, tempstring, ';');
        sidelength = stod(tempstring);
        boxSize = vector<double>({sidelength, sidelength});        
    }

        



    /***********************************************************/
    /*                  Initialize System                      */
    /***********************************************************/
    aabb::Tree simTree(2, maxDisp, period, boxSize, numParticles);
    vector<vector<double> > positions;
    vector<double> orientations;
    vector<Square> squareList;
    vector<unsigned int> tags;
    vector<unsigned int> idxs;
    vector<int> failRecord;

    if(!loadPrevious){
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
        for(int i = 0; i < positions.size(); ++i){
            file << positions[i][0] << ";" << positions[i][1] << ";" << orientations[i] << "\n";
        }
    }
    else{
        //start loading previous;
        int totalCount = numParticles * 199;
        int tempCount = 0;
        string line;
        while(tempCount < totalCount){
            getline(file, line);
            ++tempCount;
        }
        tempCount = 0;
        while(tempCount < numParticles){
            getline(file, line);
            std::istringstream sline(line);
            int col = 0;
            double previous;
            while(sline){
                std::string tempstring;
                if(!std::getline(sline, tempstring, ';')) break;
                double number = std::stod(tempstring);
                if(col == 1){
                    positions.push_back(vector<double>({previous, number}));
                    simTree.insertParticle(tempCount, positions.back(), diagRadius);
                }
                else if(col == 2){
                    orientations.push_back(number);
                    idxs.push_back(tempCount);
                    if(tempCount < activeParticles) tags.push_back(1);
                    else tags.push_back(0);
                    squareList.push_back(Square(positions.back(), number));
                }
                previous = number;
                ++col;
            }
            ++tempCount;
        }        
    }




    /***********************************************************/
    /*                  Starting Dynamics                      */
    /***********************************************************/
    
    for(int it = 0; it < nIteration; ++it){
        random_shuffle(idxs.begin(), idxs.end());
        int failed_moves = 0;
        for(unsigned int idx: idxs){
            unsigned int pt_tag = tags[idx];  //get particle id;
            double curr_orient = orientations[idx];  //get current particle orientation

            double mu = (pt_tag == 1) ? v0: 0;  //assign velocity

            double rotation_disp = myGenerator.normal_rdnm(0, DR);  //sample rotation velocity
            double disp_0 = myGenerator.normal_rdnm(mu, D0);  //sample displacement along self-propelled direction
            double disp_t = myGenerator.normal_rdnm(0, DT);   //sample transverse dispalcement 

            double disp_x = disp_0 * cos(curr_orient) - disp_t * sin(curr_orient);  //calculate dispalcement along x-axis;
            double disp_y = disp_0 * sin(curr_orient) + disp_t * cos(curr_orient);  //calculate displacement along y-axis;

            vector<double> pos({positions[idx][0] + disp_x, positions[idx][1] + disp_y});  //calculate new particle position;

            periodBoundary(pos, period, boxSize);  //adapt period boundary condition to the new particle position;

            Square tempSquare(pos, curr_orient + rotation_disp);   //create new square type;
            vector<double> lowerBound({pos[0] - diagRadius, pos[1] - diagRadius});  // lower bound of the AABB box;
            vector<double> upperBound({pos[0] + diagRadius, pos[1] + diagRadius});  // upper bound of the AABB box;
            
            aabb::AABB boundBox(lowerBound, upperBound);  //create AABB box;

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
                    }
                    else{
                        ++failed_moves;
                    }           
                }


            }
        }

        failRecord.push_back(failed_moves);

        if((it+1) % 100 == 0)
            std::cout << "Finish Iteration " << (it+1) << endl;
        if(it % stepSave == 0){
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


