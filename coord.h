#ifndef COORD
#define COORD
#include <cmath>
#include <vector>

class TField {

public:
    TField();
    void calculatedCoord(float dt, long long i);
private:
    void calculateSpeed_Up();
public:
     struct elem {
        float p[3];
    };
    struct elem speed;
    struct elem coord;
    struct elem speedUP;
    //2360591.5 == 1 //655.719861
    //5056.879 == 2 // 1.40468861
    //500 == 3 //0.138
    const float T = 0.138;
    const float g = K/(R*R);
    //6367444.7 == 1 //6367.4447
    //384400000 == 2  //384400
    const float R = 6367.4447;
    //
    const float K = 5.1625595801952e+12;

};

#endif // COORD

