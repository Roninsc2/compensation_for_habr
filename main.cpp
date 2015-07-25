#include <iostream>
#include "coord.h"
#include <fstream>

using namespace std;

ofstream fout;

TField::TField()
{
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    //speed.p[2] = 0;
}

void TField::Split(float a, int s, float& a_hi, float& a_lo)
{
    float c = (pow(2, s) + 1)*a;
    float a_big = c - a;
    a_hi = c - a_big;
    a_lo = a - a_hi;
}

float TField::TwoProduct(float a, float b, float& err)
{
    float x = a*b;
    float a_hi, a_low, b_hi, b_low;
    Split(a, 16, a_hi, a_low);
    Split(b, 16, b_hi, b_low);
    float err1, err2, err3;
    err1 = x - (a_hi*b_hi);
    err2 = err1 - (a_low*b_hi);
    err3 = err2 - (a_hi*b_low);
    err += ((a_low * b_low) - err3);
    return x;

}

void TField::CalculateSpeedUp(float err[3], float err2[3], float(& errorSpeedUp3)[3], float(& errorSpeedUp)[3])
{
    //вычисляю ускорение. vectorR - радиус вектор, вычисляется с учитываением погрешности координаты
    //vectorR = sqrt(координата_x^2 + координата_y^2) - общий вид
    float input[3] = {0,0,0};
    float errorA2[3] = {0,0,0};
    float errorA[3] = {0,0,0};
    float errorSpeedUp2[3] = {0,0,0};
    for (int l=0; l < 2; l++) {
        float temp = TwoProduct(err[l], err[l], errorA2[l]);
        temp = TwoSum(TwoProduct(err2[l], err2[l], errorA2[l]), temp, errorA[l], false);
        temp = TwoSum(temp, TwoProduct(err2[l]*2, err[l], errorA2[l]), errorA[l], false);
        temp = TwoSum(temp, TwoProduct(err2[l]*2, coord.p[l], errorA2[l]), errorA[l], false);
        temp = TwoSum(temp, TwoProduct(2*coord.p[l], err[l], errorA2[l]), errorA[l], false);
        temp = TwoSum(temp, TwoProduct(coord.p[l], coord.p[l], errorA2[l]), errorA[l], false);
        input[l] = temp;
        errorA[l] = TwoSum(errorA[l], errorA2[l], errorA2[l], true);
    }
    float temp = 0.0;
    float errorMain = 0.0;
    float inputMain =  TwoSum(errorA[0],errorA[1], errorMain, false);
    inputMain = TwoSum(input[0], inputMain, errorMain, false);
    inputMain = TwoSum(input[1], inputMain,errorMain, false);
    inputMain = TwoSum(inputMain, errorMain, errorMain, true);

    float errorVectorR = 0.0;
    float vectorR = TwoSum(sqrt(inputMain)*inputMain, (3*errorMain)*sqrt(inputMain)/2, errorVectorR, false);
    /*float x = (vectorR*errorVectorR) / (errorVectorR*errorVectorR + 2*vectorR*errorVectorR + vectorR*vectorR);
    x = -1 * K * x;*/
    if (vectorR == 0) {
        speedUP.p[0] = 0;
        speedUP.p[1] = 0;
        //speedUP.p[2] = 0;
    } else {
        //раскадываем ускорение по координатам и учитываем погрешность каждые 10 шагов
        float vectorA  = 0.0;
        float errorVectorA = 0.0;
        vectorA = TwoSum(-K/vectorR, errorVectorR*K/(vectorR*vectorR), errorVectorA, false);
        for (int l = 0; l < 2; l++) {
             errorSpeedUp2[l] = 0.0;
             errorSpeedUp[l] = 0.0;
             errorSpeedUp3[l] = 0.0;
             temp = TwoProduct(vectorA, err[l], errorSpeedUp2[l]);
             temp = TwoSum(temp, TwoProduct(vectorA, err2[l], errorSpeedUp2[l]), errorSpeedUp[l], false);
             temp = TwoSum(temp, TwoProduct(vectorA, coord.p[l], errorSpeedUp2[l]), errorSpeedUp[l], false);
             //
             temp = TwoSum(temp, TwoProduct(errorVectorA, err[l], errorSpeedUp2[l]), errorSpeedUp[l], false);
             temp = TwoSum(temp, TwoProduct(errorVectorA, err2[l], errorSpeedUp2[l]), errorSpeedUp[l], false);
             temp = TwoSum(temp, TwoProduct(errorVectorA, coord.p[l], errorSpeedUp2[l]), errorSpeedUp[l], false);
             //
             speedUP.p[l] = temp;
             //
             temp = TwoSum(errorSpeedUp[l], errorSpeedUp2[l], errorSpeedUp2[l], true);
             errorSpeedUp3[l] = TwoSum(errorSpeedUp3[l],temp, errorSpeedUp[l], true);
        }
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
    //задаем погрешности для скорости, расстояние, ускорения, погрешности от погрешности ускорения
    //в принципе, ее можно выпилить из функции, если мы ее обнуляем.
    float errorSpeedUpInCoord[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};
    float error[3] = {0,0,0};
    float errorSpeed[3] = {0,0,0};
    float errorSpeedUp[3] = {0,0,0};
    float errorSpeedUp2[3] = {0,0,0};
    float error2[3] = {0,0,0};
    float errorSpeed2[3] = {0,0,0};
    float error3[3] = {0,0,0};
    float errorSpeed3[3] = {0,0,0};

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
    float temp = 0;
    for (j = 0; j < i; j++) {
        CalculateSpeedUp(error3, error, errorSpeedUp, errorSpeedUp2);
        for (int l = 0; l < 2; l++) {
            coord.p[l] = TwoSum(error3[l], coord.p[l],  error3[l], true);
        }

        for (int l = 0; l < 2; l++) {
            errorSpeedUpInCoord[l][0] = 0.0;
            errorSpeedUpInCoord[l][1] = 0.0;
            errorSpeedUpInCoord[l][2] = 0.0;
             //coord update
            temp = TwoProduct(errorSpeedUp[l], 0.5, error2[l]);
            temp = TwoProduct(temp, time, errorSpeedUpInCoord[l][0]);
            temp = TwoSum(TwoProduct(temp, time, error2[l]), TwoProduct(errorSpeedUpInCoord[l][0],time, error2[l]), error[l], false);
            input = temp;
            temp = TwoProduct(time, errorSpeedUp2[l]*0.5, errorSpeedUpInCoord[l][2]);
            temp = TwoSum(TwoProduct(temp, time, error2[l]), TwoProduct(errorSpeedUpInCoord[l][2], time, error2[l]), error[l], false);
            input = TwoSum(temp, input, error[l], false);
            input = TwoSum(input, TwoProduct(errorSpeed[l], time, error2[l]), error[l], false);
            input = TwoSum(TwoProduct(errorSpeed3[l],time, error2[l]), input, error[l], false);
            input = TwoSum(input, TwoProduct(speed.p[l],time, error2[l]),  error[l], false);
            temp = TwoProduct(speedUP.p[l], 0.5, error2[l]);
            temp = TwoProduct(temp, time, errorSpeedUpInCoord[l][1]);
            temp = TwoSum(TwoProduct(temp, time, error2[l]), TwoProduct(errorSpeedUpInCoord[l][1],time, error2[l]), error[l], false);
            input = TwoSum(input, temp, error[l], false);
            temp = TwoSum(error[l], error2[l], error2[l], true);
            error3[l] = TwoSum(error3[l],temp, error[l], true);
            coord.p[l] = TwoSum(input, coord.p[l], error3[l], false);
            //

            speed.p[l] = TwoSum(errorSpeed3[l], speed.p[l], errorSpeed3[l], true);

             //speed update
            input = TwoSum(TwoProduct(errorSpeedUp[l],time, errorSpeed2[l]),
                           TwoProduct(speedUP.p[l] ,time, errorSpeed2[l]), errorSpeed[l], false);
            input = TwoSum(input, TwoProduct(errorSpeedUp2[l], time, errorSpeed2[l]), errorSpeed[l], false);
            temp = TwoSum(errorSpeed[l], errorSpeed2[l], errorSpeed2[l], true);
            errorSpeed3[l] = TwoSum(errorSpeed3[l],temp, errorSpeed[l], true);
            speed.p[l] = TwoSum(input, speed.p[l], errorSpeed3[l], false);

            speedUP.p[l] = TwoSum(errorSpeedUp[l], speedUP.p[l], errorSpeedUp[l], true);

            errorSpeed2[l] = 0;
            error2[l] = 0;
        }

        k++;
    }
    for (int l = 0; l < 2; l++) {
        coord.p[l] = TwoSum(error3[l], coord.p[l],error3[l], true);
    }
}

int main()
{
    fout.open("float_test322.txt");
    TField field;
    fout << field.T  << "\t" << field.R << "\t" << field.K << std::endl;
    float resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << field.K << std::endl;
    float resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << resultx << "\t" << resulty;
    long long i = 500000;
    for (;  i < 10000000000 ; i += 100000) {
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
