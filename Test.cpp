#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <complex>
#include <vector>
#include <fstream>
#include <string.h>
//#include "matplotlibcpp.h"

using namespace std;

#define cbrt(x) (x * x * x)

#define TOTAL_HORIZANTAL_POINTS 10
#define TOTAL_VERTICAL_POINTS 10
#define VERTICAL_GAP 0.05f
#define C 29979.2458f

struct Point {
    float x, y, z;

    Point(): x(0), y(0), z(0){}
    Point(float iX, float iY): x(iX), y(iY), z(0){}
    Point(float iX, float iY, float iZ): x(iX), y(iY), z(iZ){}
    Point(const Point &other): x(other.x), y(other.y), z(other.z){} 

    float operator>>(Point& other) const{
        return sqrt(pow((*this).x - other.x ,2) + pow((*this).y - other.y ,2) + pow((*this).z - other.z ,2));
    }

    Point increase(float increment) const{
        Point returnVal = (*this);
        returnVal.x += increment;
        return returnVal;
    }

    Point decrease(float decrement) const{
        Point returnVal = (*this);
        returnVal.x -= decrement;
        return returnVal;
    }
};

struct Measurement {
    Point location;
    complex<float> *s21;
    float *freq;
    size_t size;

    Measurement(): s21(nullptr), freq(nullptr){}

    Measurement(int horizantal, int vertical){

        //Read S21 parameters
        ostringstream filename;
        filename << to_string(horizantal) << "-" << to_string(vertical) << ".txt";
        readParameters(filename.str());
        readFreq("freq.txt");
        readLocation("locations.txt", horizantal, vertical);
    
    }

    void readParameters(string filename){

        ifstream in;
        in.open(filename);
        
        int count = 0;
        string line;
        while (getline(in, line))
            count++;

        size = count;
        this->s21 = new complex<float>[count];

        in.clear();
        in.seekg(0, ios::beg);

        size_t pos;
        count = 0;
        while (getline(in, line)){
            pos = line.find(" ");
            s21[count] = {stof(line.substr(0, pos)), stof(line.substr(pos+1))};
            count++;
        };
        in.close();

    }

    void readFreq(string filename){
        ifstream in;
        in.open(filename);
        string line;
        this->freq = new float[size];
        int count=0;
        while (getline(in, line)){
            freq[count] = stof(line.substr(0, line.size()-4));
            count++;
        };
        in.close();
    }

    void readLocation(string filename, int hNum, int vNum) {

        ifstream in;
        in.open(filename);
        
        size_t pos;
        int count = 0;
        string line;
        while (getline(in, line)){
            if (count == hNum){
                pos = line.find(" ");
                location = {stof(line.substr(0, pos)), stof(line.substr(pos+1)), VERTICAL_GAP*vNum};
                break;
            }
            count++;
        }

        in.close();
    }

    Measurement& operator=(Measurement &&other){
        if (s21 != nullptr){
            delete [] s21;
        }
        s21 = other.s21;
        freq = other.freq;
        size = other.size;
        location = other.location;

        other.s21 = nullptr;
        other.freq = nullptr;
        return *this;
    }

    ~Measurement(){
        if (s21 != nullptr){
            delete [] s21;
        }
        if (freq != nullptr){
            delete [] freq;
        }
        
    }
};

class TargetArea {
   public:
    uint32_t imageWidth, imageHeight, imageDepth;
    vector<Point> points;
    TargetArea(){}
    TargetArea(Point p1, Point p2, double cellSize) {
        if (cellSize == 0) throw "Cell size can not be 0.";
        if (p1.x > p2.x) swap(p1.x, p2.x);
        if (p1.y > p2.y) swap(p1.y, p2.y);
        if (p1.z > p2.z) swap(p1.z, p2.z);
        
        imageWidth = (p2.x - p1.x) / cellSize;
        imageHeight = (p2.y - p1.y) / cellSize;
        imageDepth = (p2.z - p1.z) / cellSize;

        if (p2.z == p1.z){
            points.reserve(imageWidth*imageHeight);
        } else {
            points.reserve(imageWidth*imageHeight*imageDepth);
        }
        
        for (float x = p1.x; (p2.x - x) >= -numeric_limits<float>::epsilon(); x += cellSize) {
            for (float y = p1.y; (p2.y - y) >= -numeric_limits<float>::epsilon(); y += cellSize) {
                for (float z = p1.z; (p2.z - z) >= -numeric_limits<float>::epsilon(); z += cellSize) {
                    points.push_back({x, y, z});
                    cout << x << " " << y << " " << z << endl;
                }                  
            }
        }
    }

    TargetArea& operator=(const TargetArea &other){
        points = other.points;
        return (*this);
    }

    int size(){
       return points.size();
    }
};

class Algorithm {
    public:
        virtual complex<float> * createImageMatrix(TargetArea &targetArea, const Measurement *measurements) = 0;
};

class DelayAndSum: public Algorithm {
    public:
        virtual complex<float> *createImageMatrix(TargetArea &targetArea, const Measurement *measurements) {
            complex<float> J={0, 1}, *image;
            image = new complex<float>[targetArea.size()];
            for(int k=0; k<TOTAL_HORIZANTAL_POINTS; k++){
                for(int i=0; i<measurements[k].size; i++){
                    for(int l=0; l<targetArea.size(); l++){
                        float dist1 = measurements[k].location.decrease(0.464)>>targetArea.points[l];
                        float dist2 = measurements[k].location.increase(0.464)>>targetArea.points[l];
                        image[l]+= measurements[k].s21[i]*exp(-J*measurements[k].freq[i]*2.0f*(float)M_PI*(dist1+dist2)/C);
                    }
                }
            }
            return image;
        };
};

class ImageConstructor {
    private:
        TargetArea _targetArea;
        Measurement* _measurements;
        Algorithm* _algorithm;
        complex<float>* _image;
    
    public:
        ImageConstructor(TargetArea &targetArea, int verticalPoint){
            _targetArea = targetArea;
            _measurements = new Measurement[TOTAL_HORIZANTAL_POINTS];
            for (int i = 0; i < TOTAL_HORIZANTAL_POINTS; i++){
                _measurements[i] = Measurement(i, verticalPoint);
                for(int k=0; k<_measurements[i].size; k++){
                    cout << _measurements[i].freq[k] << endl;
                }
            }
        }
        ~ImageConstructor(){
            delete [] _image;
            delete [] _measurements;
        }

        void setAlgorithm(Algorithm* algorithm) {
            _algorithm = algorithm;
        }

        void createImage() {
            _image = _algorithm->createImageMatrix(_targetArea, _measurements);
        }
};

int main(){
    TargetArea targetArea(Point(0, 2, 0), Point(2, 8, 0), 0.05f);
    DelayAndSum delayAndSum;
    ImageConstructor ImageConstructor(targetArea, 0);
    ImageConstructor.setAlgorithm(&delayAndSum);
    ImageConstructor.createImage();
    return 0;
}