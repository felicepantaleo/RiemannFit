#include <Eigen/Core>
#include <Eigen/Eigenvalues>

    #include <iostream>

    #include <iomanip>
    #include "riemann_filter.h"

    bool debug=0;
    bool debug2=0;

    #define halfpi 3.141593f/2.f
    #define c_speed  299792458
    #define max_nop 8

using namespace std;
using namespace Eigen;

namespace {

constexpr float b = 1.f;
constexpr float d = 1.e4f;
int nop;

typedef Matrix<float, Dynamic, Dynamic, 0, max_nop, max_nop> MatrixNf;
typedef Matrix<float, Dynamic, Dynamic, 0, 2*max_nop, 2*max_nop> Matrix2Nf;
typedef Matrix<float, Dynamic, Dynamic, 0, 3*max_nop, 3*max_nop> Matrix3Nf;
typedef Matrix<float, Dynamic, 1, 0, max_nop, 1> VectorNf;
typedef Matrix<float, Dynamic, 1, 0, 2*max_nop, 1> Vector2Nf;
typedef Matrix<float, Dynamic, 1, 0, 3*max_nop, 1> Vector3Nf;
typedef Matrix<float, 2, Dynamic, 0, 2, max_nop> Matrix2xNf;
typedef Matrix<float, 3, Dynamic, 0, 3, max_nop> Matrix3xNf;
typedef Matrix<float, 5, 5> Matrix5f;
typedef Matrix<float, 5, 1> Vector5f;

struct circle_fit{
Vector3f par;
Matrix3f cov;
int charge;
float chi2;
};

struct line_fit{
Vector2f par;
Matrix2f cov;
float chi2;
};

struct helix_fit{
Vector5f par;
Matrix5f cov;
int charge;
float chi2;
};

struct scatter{
    float p;
    float theta;
    float X;
};

inline float sqr(float a){
return a*a;
}

//only barrel ! TO FIX
MatrixNf Scatter_cov_rad (Matrix2xNf p2D, scatter MS){
    VectorNf rad = (p2D.row(0).array().square() + p2D.row(1).array().square()).sqrt();
    MatrixNf scatter_cov_rad = MatrixXf::Zero(nop,nop);
    float sig2 = sqr(0.015/MS.p *sqr(MS.X/cos(MS.theta))*(1+0.038*log(MS.X/cos(MS.theta))));
    for(int k=0; k<nop; k++){
        for(int l=k; l<nop; l++){
            for(int i=0; i<min(k,l); i++){
                scatter_cov_rad(k,l) += (rad(k)-rad(i))*(rad(l)-rad(i))*sig2/sqr(sin(MS.theta));
                scatter_cov_rad(l,k) = scatter_cov_rad(k,l);
            }
        }
    }
    return scatter_cov_rad;
}

Matrix2Nf cov_radtocart (Matrix2xNf p2D, MatrixNf cov_rad){
    Matrix2Nf cov_cart=MatrixXf::Zero(2*nop,2*nop);
    VectorNf rad = (p2D.row(0).array().square() + p2D.row(1).array().square()).sqrt();
    for(int k=0; k<nop; k++){
        for(int l=k; l<nop; l++){
            cov_cart(k,l) = cov_rad(k,l)*p2D(1,k)/rad(k)*p2D(1,l)/rad(l);
            cov_cart(k+nop,l+nop) = cov_rad(k,l)*p2D(0,k)/rad(k)*p2D(0,l)/rad(l);
            cov_cart(k,l+nop) = -cov_rad(k,l)*p2D(1,k)/rad(k)*p2D(0,l)/rad(l);
            cov_cart(k+nop,l) = -cov_rad(k,l)*p2D(0,k)/rad(k)*p2D(1,l)/rad(l);
            cov_cart(l,k) = cov_cart(k,l);
            cov_cart(l+nop,k+nop) = cov_cart(k+nop,l+nop);
            cov_cart(l+nop,k) = cov_cart(k,l+nop);
            cov_cart(l,k+nop) = cov_cart(k+nop,l);
        }
    }
    return cov_cart;
}

MatrixNf cov_carttorad (Matrix2xNf p2D, Matrix2Nf cov_cart){ //only diagonal terms
    MatrixNf cov_rad=MatrixXf::Zero(nop,nop);
    VectorNf rad2 = p2D.row(0).array().square() + p2D.row(1).array().square();
    for(int i=0; i<nop; i++){
        cov_rad(i,i) = (cov_cart(i,i)*sqr(p2D(1,i)) + cov_cart(i+nop,i+nop)*sqr(p2D(0,i))
            - 2.f*cov_cart(i,i+nop)*p2D(0,i)*p2D(1,i))/rad2(i);
        if(rad2(i)<1e-4) cov_rad(i,i) = cov_cart(i,i); //TO FIX
    }
    return cov_rad;
}

VectorNf Weight(Matrix2xNf p2D, MatrixNf cov_rad_inv){
    VectorNf weight(nop);
    for(int i=0; i<nop; i++){
        weight(i) = cov_rad_inv.col(i).sum();
    }
    return weight;
}

float chi2_circle(Matrix2xNf p2D, Matrix2Nf V, Vector3f par_uvr){
    float chi2 = 0;
    for(int i=0; i<nop; i++){
        float x_ = p2D(0,i)-par_uvr(0);
        float y_ = p2D(1,i)-par_uvr(1);
        float x_2 = sqr(x_);
        float y_2 = sqr(y_);
        chi2 += sqr(sqrt(x_2+y_2)-par_uvr(2))/((V(i,i)*x_2+V(i+nop,i+nop)*y_2+2*V(i,i+nop)*x_*y_)/(x_2+y_2));
    }
    return chi2;
}

float chi2_line(Matrix2xNf p2D, VectorNf y_err2, Vector2f par_line){
    float chi2 = 0;
        for(int i=0; i<nop; i++){
        chi2 += sqr(p2D(1,i)-p2D(0,i)*par_line(0)-par_line(1)) / (1+sqr(par_line(0)))
                / (y_err2(i)*sqr(cos(atan(par_line(0)))) + y_err2(i)*sqr(sin(atan(par_line(0)))));
    }
    return chi2;
}

inline int Charge (Matrix2xNf p2D, Vector3f par_uvr) { //error to be computed TO FIX
    float dir = (p2D(0,1)-p2D(0,0))*(par_uvr(1)-p2D(1,0))-(p2D(1,1)-p2D(1,0))*(par_uvr(0)-p2D(0,0));
    return (dir > 0) ? -1 : 1;
}

Vector3f par_transformation(Vector3f par_uvr, int charge, float B_field){
    Vector3f par_pak;
    float phi = (charge > 0) ? atan2(par_uvr(0), -par_uvr(1)) : atan2(-par_uvr(0), par_uvr(1));
    par_pak <<  phi,
                charge * (sqrt(sqr(par_uvr(0))+sqr(par_uvr(1)))-par_uvr(2)),
                par_uvr(2)*B_field;
    return par_pak;
}

// return the eigenvector associated to the minimum eigenvalue
Vector3f min_eigen3D(Matrix3f A){
    EigenSolver<Matrix3f> solver(A); // evaluate eigenvalues and eigenvector
    Vector3cf lambdac = solver.eigenvalues(); //why can't I cast here ??
    Matrix3cf eigenvectc = solver.eigenvectors();
    Vector3f lambda = lambdac.real().cast<float>(); //check if real ! TO FIX
    Matrix3f eigenvect = eigenvectc.real().cast<float>();
    int minindex =0;
    lambda.minCoeff(&minindex);
    Vector3f n = eigenvect.col(minindex);
    return n;
}

Vector2f min_eigen2D(Matrix2f A){
    EigenSolver<Matrix2f> solver(A); // evaluate eigenvalues and eigenvector
    Vector2cf lambdac = solver.eigenvalues(); //why can't I cast here ??
    Matrix2cf eigenvectc = solver.eigenvectors();
    Vector2f lambda = lambdac.real().cast<float>(); //check if real ! TO FIX
    Matrix2f eigenvect = eigenvectc.real().cast<float>();
    int minindex =0;
    lambda.minCoeff(&minindex);
    Vector2f n = eigenvect.col(minindex);
    return n;
}

//     a       ||   0|   1|   2|   3|   4|   5
// nu(a)=(i,j) || 0,0| 0,1| 0,2| 1,1| 1,2| 2,2
inline int nu(int a, int* i, int* j){
    *i = (a==0 || a==1 || a==2) ? 0 : (a==3 || a==4) ? 1 : (a==5) ? 2 : 3;
    *j = (a==2 || a==4 || a==5) ? 2 : (a==1 || a==3) ? 1 : (a==0) ? 0 : 3;
    if (*i == 3) return 1;
    return 0;
}



circle_fit Circle_fit(Matrix2xNf p2D, Matrix2Nf V, float B_field,
                      bool return_err = true, bool scattering=true, scatter* MS=nullptr){

    //INITIALIZATION

    //SORTING !! TO FIX
    Matrix2xNf p2D_ = p2D;
    Matrix3xNf p3D(3,nop);

    MatrixNf cov_rad = cov_carttorad(p2D_, V);
    if(scattering && MS != nullptr){
        MatrixNf scatter_cov_rad = Scatter_cov_rad(p2D_, *MS);
        V += cov_radtocart(p2D_, scatter_cov_rad);
        cov_rad += scatter_cov_rad;
    }
    MatrixNf G = cov_rad.inverse();
    G /= G.sum();
    VectorNf weight = Weight(p2D_, G);
    if(debug){
        cout << "cov_rad:\n" << cov_rad << endl << endl;
        cout << "G:\n" << G << endl << endl;
        cout << "weight:\n" << weight.transpose() << endl;
    }


    //CENTER & SCALE 2D POINTS
    float umean = p2D.row(0).mean();
    float vmean = p2D.row(1).mean();
    p2D.row(0) = p2D.row(0).array() - umean;
    p2D.row(1) = p2D.row(1).array() - vmean;
    Vector2Nf mc(2*nop);
    mc << p2D.row(0).transpose(), p2D.row(1).transpose(); //useful for error propagation
    float q = p2D.array().square().sum();
    float s = b*sqrt(nop/q); //scaling factor (b is an arbitrary constant)
    p2D *= s;

    //CALCULATE TRASFORMED POINTS IN 3D
    p3D.block(0,0,2,nop) = p2D;
    p3D.row(2) = p2D.row(0).array().square() + p2D.row(1).array().square();

    //CALCULATE & MINIMIZE COST FUNCTION
    Matrix3f A = Matrix3f::Zero();
    Vector3f r0 = p3D * weight;

    Matrix3xNf temp = p3D - r0*RowVectorXf::Constant(nop,1.f);
    A = temp*G*temp.transpose();

    Vector3f n = min_eigen3D(A);
    n *= (n(2)>0) ? 1 : -1;
    float c = -n.transpose()*r0;

    //CALCULATE CIRCUMFERENCE PARAMETER
    Vector3f par_uvr_;
        par_uvr_ << -n(0)/(2.f*n(2)),
                    -n(1)/(2.f*n(2)),
                    sqrt((1.f-sqr(n(2))-4.f*c*n(2))/(4.f*sqr(n(2))));
    Vector3f par_uvr;
        par_uvr <<  par_uvr_(0)/s + umean, par_uvr_(1)/s + vmean, par_uvr_(2)/s;

    circle_fit circle;
    circle.charge = Charge(p2D_, par_uvr);
    circle.par = par_transformation(par_uvr, circle.charge, B_field);
    if(!scattering || (scattering && MS!=nullptr))
        circle.chi2 = chi2_circle(p2D_, V, par_uvr);



    return circle;
}

line_fit Line_fit (Matrix3xNf p3D, Matrix3Nf V, circle_fit circle,
                   float B_field, bool return_err=true){

    //INITIALIZATION
    Matrix2xNf p2D(2,nop);

    //VectorNf x_err2 = X_err(p3D, V, circle);
    VectorNf y_err2 = V.block(2*nop,2*nop,nop,nop).diagonal();
    //float k = x_err2.array().sqrt().mean()/y_err2.array().sqrt().mean();
    //if(debug) cout << "k:  " <<  k <<endl;

    VectorNf weight = 1.f/(y_err2).array();
    weight /= weight.sum();

    //CALCULATE TRASFORMED POINTS IN 2D
    p2D.row(1) = p3D.row(2);
    //p2D.row(0) = p3D.row(1);
    Matrix<float, 1, Dynamic> ciccio=((p3D.row(0).array() - cos(circle.par(0)-halfpi)*circle.par(1)).square()
        + (p3D.row(1).array() - sin(circle.par(0)-halfpi)*circle.par(1)).square()).sqrt()
        / (circle.par(2)/B_field*2.f);
        for(int i=0; i<nop; i++){
            if(ciccio(i) < -1) ciccio(i)=-1.f;
            else if(ciccio(i) >1) ciccio(i)=1.f;
        }

    p2D.row(0) = 2.f*asin(ciccio.array())*(circle.par(2)/B_field);

    //CALCULATE & MINIMIZE COST FUNCTION
    Matrix2f A = Matrix2f::Zero();
    Vector2f r0 = p2D * weight;

    for (int i=0; i<nop; i++) A += weight(i)*((p2D.col(i)-r0)*(p2D.col(i)-r0).transpose());

    Vector2f n = min_eigen2D(A);
    float c = -n.transpose()*r0;

    //CALCULATE LINE PARAMETER
    line_fit line;
    line.par << -n(0)/n(1), -c*sqrt(sqr(n(0))+sqr(n(1)))/n(1);
    line.chi2 = chi2_line(p2D, y_err2, line.par);

    //ERROR PROPAGATION
    if(return_err){
        //auxiliary quantities
        float sig2 = y_err2.mean(); //TO FIX
        float S = (A(0,0) + A(1,1))*nop;
        float n0_2 = sqr(n(0));
        float n1_2 = sqr(n(1));
        float sqrt_ = sqrt(n1_2+n0_2);
        float x_ =  p2D.row(0).sum()/nop;
        float y_ =  p2D.row(1).sum()/nop;
        float corr = sqr(1.131); //TO FIX
        float C13 = sig2*n(1)*(n(0)*y_-n(1)*x_)/S;
        float C23 =-sig2*n(0)*(n(0)*y_-n(1)*x_)/S;
        float C33 = corr*sig2*sqr(1/nop + n(0)*y_-n(1)*x_)/S;
        Matrix3f C;
        C <<    sig2*n1_2/S, -sig2*n(0)*n(1)/S, C13,
                -sig2*n(0)*n(1)/S, sig2*n0_2/S, C23,
                C13, C23, C33;

        Matrix<float, 2, 3> J;
        J <<    -1.f/n(1), n(0)/n1_2, 0,
                -c*n(0)/(n(1)*sqrt_), n0_2*c/(n1_2*sqrt_), -sqrt_/n(1);

        line.cov = J * C * J.transpose();
    }

    return line;
}

helix_fit Helix_fit(Matrix3xNf hits, Matrix3Nf hits_cov, float B_field,
                    bool return_err=true, bool scattering=true){
    nop = hits.cols();
    circle_fit circle;
    line_fit line;
    scatter MS;
    if(scattering){
        circle = Circle_fit(hits.block(0,0,2,nop), hits_cov.block(0,0,2*nop,2*nop),
                            B_field, false, true, nullptr);
        line = Line_fit(hits, hits_cov, circle, B_field, false);
        MS.theta = atan(1/line.par(0));
        MS.p = circle.par(2)/cos(MS.theta);
        MS.X = 0.04f;
    }
    circle = Circle_fit(hits.block(0,0,2,nop), hits_cov.block(0,0,2*nop,2*nop),
                        B_field, return_err, scattering, &MS);
    line = Line_fit(hits, hits_cov, circle, B_field, return_err);
    helix_fit helix;
    helix.par << circle.par, line.par;
    helix.cov = MatrixXf::Zero(5,5);
    if(return_err){
        helix.cov.block(0,0,3,3) = circle.cov;
        helix.cov.block(3,3,2,2) = line.cov;
    }
    helix.charge = circle.charge;
    helix.chi2 = circle.chi2 + line.chi2;
    return helix;
}

}



bool havesamephi(Vector3f p1, Vector3f p2, float phiCut){
    Vector2f p1_ = p1.block(0,0,2,1);
    Vector2f p2_ = p2.block(0,0,2,1);
    return (((p1_.transpose()*p2_)/(p1_.norm()*p2_.norm())).norm() > phiCut);
}

bool areAlignedRZ(const float r1, const float z1,const float r2,const float z2,const float r3,const float z3, float ptmin, float thetaCut)
    {
        float radius_diff = std::abs(r1 - r3);
        float distance_13_squared = radius_diff*radius_diff + (z1 - z3)*(z1 - z3);

        float pMin = ptmin*std::sqrt(distance_13_squared); //this needs to be divided by radius_diff later

        float tan_12_13_half_mul_distance_13_squared = fabs(z1 * (r2 - r3) + z2 * (r3 - r1) + z3 * (r1 - r2)) ;
        return tan_12_13_half_mul_distance_13_squared * pMin <= thetaCut * distance_13_squared * radius_diff;
    }

void findTps(const vector< event_t>& events, vector < vector <vector <array <int,4>>>>& tracking_particles){
    for (int event_id=0; event_id <events.size(); event_id++){
        for(unsigned int offset=0; offset<4; offset++){
            int first_lid = offset;
            for(unsigned int i=0; i<events[event_id][first_lid].size(); i++){
                if(events[event_id][first_lid][i].inner_tp == events[event_id][first_lid][i].outer_tp){
                    for(unsigned int j=0; j<events[event_id][first_lid+1].size(); j++){
                        if( events[event_id][first_lid][i].outer_tp == events[event_id][first_lid+1][j].inner_tp &&
                            events[event_id][first_lid+1][j].inner_tp == events[event_id][first_lid+1][j].outer_tp){
                            for(unsigned int k=0; k<events[event_id][first_lid+2].size(); k++){
                                if(events[event_id][first_lid+1][j].outer_tp == events[event_id][first_lid+2][k].inner_tp &&
                                    events[event_id][first_lid+2][k].inner_tp == events[event_id][first_lid+2][k].outer_tp){
                                        tracking_particles[event_id][first_lid].emplace_back(array<int, 4>{{events[event_id][first_lid][i].inner_tp, i, j, k}});
                                    //cout << tracking_particles[event_id][first_lid].back()[0] << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void track_part(const vector <event_t> &events,  vector <vector <vector <array<int,4>>>> &tracking_particles, ofstream & myfile){
    findTps(events, tracking_particles);
    for (int ev_id=0; ev_id<tracking_particles.size(); ev_id++){

        for(unsigned int offset=0; offset<4; offset++){
            int first_lid = offset;
             sort(tracking_particles[ev_id][offset].begin(),tracking_particles[ev_id][offset].end(),[] (const array<int,4>& a, const array<int,4>& b ) { return a[0]< b[0]; } );
             tracking_particles[ev_id][offset].erase(unique(tracking_particles[ev_id][offset].begin(),tracking_particles[ev_id][offset].end(),
                                                            [] (const array<int,4>& a, const array<int,4>& b ) { return a[0]== b[0]; } ),
                                                    tracking_particles[ev_id][offset].end());
                    }
                }

    for (int ev_id=0; ev_id<tracking_particles.size(); ev_id++){
        for(unsigned int offset=0; offset<4; offset++){
            int first_lid = offset;
            for(unsigned int tp=0; tp<tracking_particles[ev_id][offset].size(); tp++){
             myfile << ev_id << "  " << first_lid << " " << tracking_particles[ev_id][offset][tp][0] << endl;
            }
        }
    }
}

int main(){
    ofstream myfile("tracking_particles.txt", ofstream::out);
    ofstream comecazzovuoi("fileditesto.txt", ofstream::out);
    vector <event_t> events;
    vector <vector <vector <array<int,4>>>> tracking_particles;
    ciccio("data.csv", events);
    Matrix<float, 3, 7> hit;
    Matrix<float, 21, 21> cov=MatrixXf::Zero(21,21);
    Vector3f zero; zero << 0,0,0;
    unsigned long long int counter=0;
    tracking_particles.resize(events.size());
    for(auto& tp:tracking_particles){
        tp.resize(4);
    }

    track_part(events, tracking_particles, myfile);

 for (auto& ev: events){

        for(unsigned int offset=0; offset<4; offset++){
            int first_lid = offset;
            //int first_lid=0;
            //int i=3194; int j=1363; int k=576;
            for(unsigned int i=0; i<ev[first_lid].size(); i++){
                if(i%10==0){
                unsigned long long int index = (unsigned long long int)i*ev[first_lid+1].size()*ev[first_lid+2].size();
                cout << "processed: " << counter << "    " <<  index << endl;
                }
                float r1 = sqrt(sqr(ev[first_lid][i].point1(0))+sqr(ev[first_lid][i].point1(1)));
                for(unsigned int j=0; j<ev[first_lid+1].size(); j++){
                    if(havesamephi(ev[first_lid][i].point1,ev[first_lid+1][j].point1,0.9f)){
                        float r2 = sqrt(sqr(ev[first_lid+1][j].point1(0))+sqr(ev[first_lid+1][j].point1(1)));

                        for(unsigned int k=0; k<ev[first_lid+2].size(); k++){
                            if(havesamephi(ev[first_lid+1][j].point1,ev[first_lid+2][k].point1,0.9f)){
                                float r3 = sqrt(sqr(ev[first_lid+2][k].point1(0))+sqr(ev[first_lid+2][k].point1(1)));

                                if(areAlignedRZ(r1, ev[first_lid][i].point1(2), r2, ev[first_lid+1][j].point1(2),r3, ev[first_lid+2][k].point1(2),0.5f, 0.01)){

                                    counter++;
                                    hit <<  zero, ev[first_lid][i].point1, ev[first_lid][i].point2, ev[first_lid+1][j].point1,
                                            ev[first_lid+1][j].point2,ev[first_lid+2][k].point1, ev[first_lid+2][k].point2;

                                    cov(0,0)=0.1;
                                    cov(1,1)=0.1;
                                    cov(2,2)=5;
                                    for(int z=0; z < 3; z++){
                                         cov(6*z+3,6*z+3) = ev[first_lid][i].err_point1(z);
                                         cov(6*z+4,6*z+4) = ev[first_lid][i].err_point2(z);
                                         cov(6*z+5,6*z+5) = ev[first_lid+1][j].err_point1(z);
                                         cov(6*z+6,6*z+6) = ev[first_lid+1][j].err_point2(z);
                                         cov(6*z+7,6*z+7) = ev[first_lid+2][k].err_point1(z);
                                         cov(6*z+8,6*z+8) = ev[first_lid+2][k].err_point2(z);
                                    }
                                    cov.array().square();

                                    helix_fit helix = Helix_fit(hit, cov, 3.8*c_speed/pow(10,9)/100, false, false);
                                    //cout << i << " " << j << "  " << k << endl;
                                    if(helix.chi2 < 100){
                                        //cout << ev[first_lid][i].inner_tp   << "   " << helix.chi2 << endl;
                                        comecazzovuoi  << ev[first_lid][i].inner_tp   << "   "
                                                << ev[first_lid][i].outer_tp    << "   "
                                                << ev[first_lid+1][j].inner_tp  << "   "
                                                << ev[first_lid+1][j].outer_tp  << "   "
                                                << ev[first_lid+2][k].inner_tp  << "   "
                                                << ev[first_lid+2][k].outer_tp  << "   " << helix.chi2 << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

comecazzovuoi.close();
myfile.close();


    return 0;
}
