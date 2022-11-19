/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#ifndef LA_H
#define LA_H

#include <cmath>
#include <vector>

using namespace std;
typedef vector<float> V1;
typedef vector<vector<float> > V2;

V1  rotMatrix(V1 vec, V1 axis, float theta);
V1 const_divide(float par,V1 c);
float _distance(V1 &p1,V1 &p2);
void _multiply(V1 &c1,V1 &c2);
V1 _add(V1 c1, V1 c2);
void _subtract(V1 &c1, V1 &c2, V1 &c3);
V1 const_multiply(float par,V1 c);
void _normalize(V1 &c);
void cross_product(V1 &c1, V1 &c2, V1 &c3);
float dot_product(V1 &c1,V1 &c2);
float calcTorsionAngle(V1 &c1, V1 &c2, V1 &c3, V1 &c4);
V1 class1(V1 &a, V1 &b, V1 &c, float bond_len, float di_angle, float theta);
V1 calcCoordinate(V1 &a, V1 &b, V1 &c, float bond_len, float theta, float di_angle);
V1 recalcCoordinate(V1 b, V1 c, V1 d, float bond_len);
V2 class2( V1 pos_CA, V1 pos_CG, V1 pos_CB, float bond_len = 0.0);
V1 class3(V1 pos_CB, V1 pos_N, V1 pos_CA);
V1 class4(V1 pos_NE, V1 pos_CZ, V1 pos_NH1);
V1 class5(V1 pos_CD, V1 pos_NE, V1 pos_CZ, float bond_len);
V1 class6_Ser(V1 pos_OG, V1 pos_CB, V1 HB1);
V1 class6_Thr(V1 pos_OG1, V1 pos_CB, V1 pos_CG2);
V1 class6_Tyr(V1 pos_OH, V1 pos_CZ, V1 pos_CE2);
V1 class6_Cys(V1 pos_SG, V1 pos_CB, V1 pos_CA); 

#endif