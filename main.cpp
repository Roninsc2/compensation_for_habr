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
}

void TField::calculateSpeed_Up() {

    float vectorR = sqrt(coord.p[0]*coord.p[0] + coord.p[1]*coord.p[1]);
    float vectorA = -K/(vectorR * vectorR * vectorR);
    for (int l = 0; l < 2; l++) {
         speedUP.p[l] = vectorA * coord.p[l];
    }
}

void TField::calculatedCoord(float time, long long i) {
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    long long j = 0;
    speedUP.p[0] = 0;
    speedUP.p[1] =0;
    speedUP.p[2] = 0;
    for (j = 0; j < i; j++) {
        calculateSpeed_Up();
        for (int l = 0; l < 2; l++) {
            coord.p[l] += speed.p[l] * time + 0.5*speedUP.p[l]*time*time;
            speed.p[l] += speedUP.p[l] * time;
        }
    }
}

int main() {
    fout.open("float_test32_simple.txt");
    TField field;
    fout << field.T  << "\t" << field.R << "\t" << field.K << std::endl;
    float resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cout << field.K << std::endl;
    float resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cout << resultx << "\t" << resulty;
    long long i = 1000;
    for (;  i < 200000000; i+= 1000) {
        float time = field.T/(float)i;
        field.calculatedCoord(time, i);
        //fout << field.coord.p[0] << "\t" << field.coord.p[1] << "\t"; - вывод координат

        //fout  << (field.coord.p[0] - resultx)/resultx << "\t" << (field.coord.p[1] - resulty)/resulty << ";" << i << std::endl;
        //вывод относительной ошибки для каждой из координат

        fout << sqrt((-resultx + field.coord.p[0])*(-resultx + field.coord.p[0])
                  + (-resulty + field.coord.p[1])*(-resulty + field.coord.p[1]))/field.R ;
        //вывод полной относительной ошибки
        fout << ";"<< i << std::endl;

    }
    fout.close();
    return 0;
}
