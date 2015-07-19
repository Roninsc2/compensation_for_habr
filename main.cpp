#include <iostream>
#include "coord.h"
#include <fstream>

using namespace std;

ofstream fout;

TField::TField() {
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    //speed.p[2] = 0;
}

void TField::CalculateSpeedUp(float err[3], float(& errorSpeedUp)[3]) {
    float vectorR = sqrt(coord.p[0]*coord.p[0] + coord.p[1]*coord.p[1] + (err[0]*err[0] + err[1]*err[1]));
    float vectorA = -K/(vectorR * vectorR * vectorR);
    for (int l = 0; l < 2; l++) {
         errorSpeedUp[l] = 0.0;
         speedUP.p[l] = TwoSum(vectorA * coord.p[l], err[l] * vectorA, errorSpeedUp[l], false);
    }
}
//компенсация суммирования по Шевчуку
float TField::TwoSum(float a, float b, float& error,  bool isNull) {
    float x = a + b;
    float b_virt = x - a;
    float a_virt = x - b_virt;
    float b_roundoff = b - b_virt;
    float a_roudnoff = a - a_virt;
    float y = a_roudnoff + b_roundoff;
    if (!isNull) {
        error += y;
    } else {
        error = y;
    }
    return x;
}


void  TField::CalculatedCoord(float time, long long i) {
    float error[3] = {0,0,0};
    float errorSpeedUp[3] = {0,0,0};
    float errorSpeed[3] = {0,0,0};

    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    float input = 0;
    long long j = 0;
    speedUP.p[0] = 0;
    speedUP.p[1] =0;
    speedUP.p[2] = 0;
    int k = 0;
    for (j = 0; j < i; j++) {
        CalculateSpeedUp(error, errorSpeedUp);
        if (k == 10) {
            for (int l = 0; l < 2; l++) {
                coord.p[l] = TwoSum(coord.p[l], error[l], error[l], true);
                speed.p[l] = TwoSum(speed.p[l], errorSpeed[l], errorSpeed[l], true);
            }
            k = 0;
        }
        for (int l = 0; l < 2; l++) {
            input = errorSpeed[l]*time + 0.5*errorSpeedUp[l]*time*time + speed.p[l]*time + speedUP.p[l]*0.5*time*time;
            coord.p[l] = TwoSum(input, coord.p[l], error[l], false);
            input = errorSpeedUp[l]*time + speedUP.p[l]*time;
            speed.p[l] = TwoSum(input, speed.p[l], errorSpeed[l], false);
        }
        k++;

    }
    for (int l = 0; l < 2; l++) {
        coord.p[l] = TwoSum(error[l], coord.p[l],error[l], true);
    }
}

int main() {
    fout.open("float_test31.txt");
    TField field;
    fout << field.T  << "\t" << field.R << "\t" << field.K << std::endl;
    float resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << field.K << std::endl;
    float resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << resultx << "\t" << resulty;
    long long i = 100000;
    for (;  i < 10000000000; i+= 1000) {
        float time = field.T/(float)i;
        field.CalculatedCoord(time, i);
        //fout << field.coord.p[0] << "\t" << field.coord.p[1] << "\t";
        //fout  << (field.coord.p[0] - resultx)/resultx << "\t" << (field.coord.p[1] - resulty)/resulty << ";" << i << std::endl;
        fout << sqrt((-resultx + field.coord.p[0])*(-resultx + field.coord.p[0])
                  + (-resulty + field.coord.p[1])*(-resulty + field.coord.p[1]))/field.R ;
        fout << ";"<< i << std::endl;
    }
    fout.close();
    return 0;
}
