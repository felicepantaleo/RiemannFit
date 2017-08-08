 #include <sstream>
 #include <fstream>
 #include <Eigen/Core>
 #include <vector>
 #include <string>


using namespace Eigen;

struct data{
Vector3f point1;
Vector3f err_point1;
Vector3f point2;
Vector3f err_point2;
int inner_tp;
int outer_tp;
float pt;
};

typedef std::vector< std::vector <data > > event_t;


 void ciccio(std::string input_file, std::vector <event_t>& events){
    char line[1024];
    std::ifstream myfile(input_file.c_str(), std::ifstream::in);
    int last_event_id = -1;
    while(myfile.good() ){
        int event_id;
        int layer_id;
        myfile.getline(line, 1023);
        std::string stringline(line);
        if (stringline.empty())
            continue;
        std::stringstream ss;
        ss.str(stringline);
        data hit;
        std::string useless;
        ss >> event_id >> layer_id >> hit.point1(0) >> hit.point1(1) >> hit.point1(2) >>
        hit.err_point1(0) >> hit.err_point1(1) >> hit.err_point1(2)
        >> hit.inner_tp
        >> hit.point2(0) >> hit.point2(1) >> hit.point2(2) >>
        hit.err_point2(0) >> hit.err_point2(1) >> hit.err_point2(2) >>
        hit.outer_tp >> hit.pt;
        bool isnewevent = (last_event_id!=event_id);
        if(isnewevent){
            last_event_id=event_id;
            event_t temp_event;
            temp_event.resize(6);
            temp_event[layer_id-1].push_back(hit);
            events.push_back(temp_event);

        }
        else{
            events.back()[layer_id-1].push_back(hit);
        }
    }

 }
