/*******************************************************************************************************************************
Copyright (c) 2022 Osita Sunday Nnyigide (osita@protein-science.com)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#include "LinearAlgebra.hpp"

V1  rotMatrix(V1 vec, V1 axis, float theta)
      {
        V2 matrix = {{},{},{}};
        matrix[0] = {0.0, 0.0, 0.0};
        matrix[1] = {0.0, 0.0, 0.0};
        matrix[2] = {0.0, 0.0, 0.0};

        float axis_length = sqrt((axis[0]*axis[0] +
              axis[1]*axis[1] + axis[2]*axis[2]));

        float xNorm = axis[0]/axis_length;
        float yNorm = axis[1]/axis_length;
        float zNorm = axis[2]/axis_length;

        float sin_theta         = sin(theta);
        float cos_theta         = cos(theta);
        float one_costheta = 1.0 - cos_theta;

        matrix[0][0] = cos_theta       + xNorm*xNorm*one_costheta;
        matrix[0][1] = xNorm*yNorm*one_costheta - zNorm*sin_theta;
        matrix[0][2] = xNorm*zNorm*one_costheta + yNorm*sin_theta;
        matrix[1][0] = xNorm*yNorm*one_costheta + zNorm*sin_theta;
        matrix[1][1] = cos_theta        + yNorm*yNorm*one_costheta;
        matrix[1][2] = yNorm*zNorm*one_costheta - xNorm*sin_theta;
        matrix[2][0] = xNorm*zNorm*one_costheta - yNorm*sin_theta;
        matrix[2][1] = yNorm*zNorm*one_costheta + xNorm*sin_theta;
        matrix[2][2] = cos_theta        + zNorm*zNorm*one_costheta;

        V1 arbrot    = { vec[0]*matrix[0][0]+ vec[1]*matrix[1][0]+ vec[2]*matrix[2][0],
                         vec[0]*matrix[0][1]+ vec[1]*matrix[1][1]+ vec[2]*matrix[2][1],
                         vec[0]*matrix[0][2]+ vec[1]*matrix[1][2]+ vec[2]*matrix[2][2]
                       };

          return arbrot;
      }

V1 const_divide(float par,V1 c)
      {
        c[0]/=par;
        c[1]/=par;
        c[2]/=par;
        return c;
      }

float _distance(V1 &p1,V1 &p2)
      {
        return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+
               (p1[1]-p2[1])*(p1[1]-p2[1])+
               (p1[2]-p2[2])*(p1[2]-p2[2]));
      }

void _multiply(V1 &c1,V1 &c2)
      {
        V1 c3  = {};
        c3.clear();
        c3.push_back(c1[0]*c2[0]);
        c3.push_back(c1[1]*c2[1]);
        c3.push_back(c1[2]*c2[2]);
      }

 V1 _add(V1 c1, V1 c2)
      {
        V1 c3  = {};
        c3.push_back(c1[0]+c2[0]);
        c3.push_back(c1[1]+c2[1]);
        c3.push_back(c1[2]+c2[2]);
        return c3;
      }

 void _subtract(V1 &c1, V1 &c2, V1 &c3)
      {
        c3.clear();
        c3.push_back(c1[0]-c2[0]);
        c3.push_back(c1[1]-c2[1]);
        c3.push_back(c1[2]-c2[2]);
      }

V1 const_multiply(float par,V1 c)
      {
        c[0]*=par;
        c[1]*=par;
        c[2]*=par;
        return c;
      }

void _normalize(V1 &c)
      {
        float len=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
        c[0]/=len;
        c[1]/=len;
        c[2]/=len;
      }

void cross_product(V1 &c1, V1 &c2, V1 &c3)
      {
        c3.clear();
        c3.push_back(c1[1]*c2[2]-c1[2]*c2[1]);
        c3.push_back(c1[2]*c2[0]-c1[0]*c2[2]);
        c3.push_back(c1[0]*c2[1]-c2[0]*c1[1]);
      }

float dot_product(V1 &c1,V1 &c2)
      {
        return c1[0]*c2[0]+c1[1]*c2[1]+c1[2]*c2[2];
      }

float calcTorsionAngle(V1 &c1, V1 &c2, V1 &c3, V1 &c4)
      {

        V1 bv12, bv32, bv43, pv13, pv24, vx;

        _subtract(c1, c2, bv12);
        _subtract(c3, c2, bv32);
        _subtract(c4, c3, bv43);

        cross_product(bv12, bv32, pv13);
        cross_product(bv43, bv32, pv24);

        float projec = dot_product(pv13, pv24);
        float sqrd13 = dot_product(pv13, pv13);
        float sqrd24 = dot_product(pv24, pv24);

        float cosine  = projec / sqrt(sqrd13*sqrd24);

        if (cosine < -1.0)
          {
            cosine = -1.0;
          }

        else if (cosine > 1.0)
          {
            cosine = 1.0;
          }
        float radian  = acos(cosine);

        cross_product(pv24, bv32, vx);
        if (dot_product(pv13, vx) < 0)
          {
              radian = -radian;
          }

        return radian*180.0/3.142;
      }

V1 class1(V1 &a, V1 &b, V1 &c, float bond_len, float di_angle, float theta)
      {

        V1 u, x, v, w, q, e, u1, y1, z, n, pos, c4, f1, f2, pos_temp, pos_BL;

        di_angle = di_angle*3.142/180.0;

        _subtract(c, b, u);
        _subtract(a, b, x);

        f1 = const_multiply(dot_product(x, u)/dot_product(u,u),u);
        _subtract(x, f1, v);

        cross_product(u, x, w);

        _normalize(v);
        _normalize(w);
        q = const_multiply(cos(di_angle), v);
        e = const_multiply(sin(di_angle), w);

        pos_temp = _add(b, _add(q, e));

        _subtract(b, c, u1);
        _subtract(pos_temp, c, y1);

        float mag_y1  = sqrt((y1[0]*y1[0])+(y1[1]*y1[1])+(y1[2]*y1[2]));
        float mag_u1  = sqrt((u1[0]*u1[0])+(u1[1]*u1[1])+(u1[2]*u1[2]));

        float theta_bcd = acos(dot_product(u1,y1)/(mag_u1*mag_y1))*180/3.142;
        float rotate    = theta - theta_bcd;

        cross_product(u1, y1, z);
        _normalize(z);

        cross_product(z, y1, f2);
        pos = _add(
                    _add(c,const_multiply(cos(rotate*3.142/180), y1)),
                    _add(const_multiply(sin(rotate*3.142/180), f2),
                        (const_multiply(1.0-cos(rotate*3.142/180),
                        const_multiply(dot_product(z, y1), z)))));

        _subtract(pos, c, c4);
        pos_BL = _add(const_multiply(bond_len/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), c);

        return pos_BL;
      }

V1 calcCoordinate(V1 &a, V1 &b, V1 &c, float bond_len, float theta, float di_angle)
      {

        V1 u, x, v, w, q, e, u1, y1, z, n, pos, c4, f1, f2, pos_temp, pos_BL;

        di_angle = di_angle*3.142/180.0;

        _subtract(c, b, u);
        _subtract(a, b, x);

        f1 = const_multiply(dot_product(x, u)/dot_product(u,u),u);
        _subtract(x, f1, v);

        cross_product(u, x, w);

        _normalize(v);
        _normalize(w);
        q = const_multiply(cos(di_angle), v);
        e = const_multiply(sin(di_angle), w);

        pos_temp = _add(b, _add(q, e));

        _subtract(b, c, u1);
        _subtract(pos_temp, c, y1);

        float mag_y1  = sqrt((y1[0]*y1[0])+(y1[1]*y1[1])+(y1[2]*y1[2]));
        float mag_u1  = sqrt((u1[0]*u1[0])+(u1[1]*u1[1])+(u1[2]*u1[2]));

        float theta_bcd = acos(dot_product(u1,y1)/(mag_u1*mag_y1))*180/3.142;
        float rotate    = theta - theta_bcd;

        cross_product(u1, y1, z);
        _normalize(z);

        cross_product(z, y1, f2);
        pos = _add(
                    _add(c,const_multiply(cos(rotate*3.142/180), y1)),
                    _add(const_multiply(sin(rotate*3.142/180), f2),
                        (const_multiply(1.0-cos(rotate*3.142/180),
                        const_multiply(dot_product(z, y1), z)))));

        _subtract(pos, c, c4);
        pos_BL = _add(const_multiply(bond_len/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), c);

        return pos_BL;
      }

V1 recalcCoordinate(V1 b, V1 c, V1 d, float bond_len)
      {
        V1 u1, y1, z, pos, c4, f1, f2, pos_BL;
        float theta = 107.0;

        _subtract(b, c, u1);
        _subtract(d, c, y1);

        float mag_y1  = sqrt((y1[0]*y1[0])+(y1[1]*y1[1])+(y1[2]*y1[2]));
        float mag_u1  = sqrt((u1[0]*u1[0])+(u1[1]*u1[1])+(u1[2]*u1[2]));

        float theta_bcd = acos(dot_product(u1,y1)/(mag_u1*mag_y1))*180/3.142;
        float rotate    = theta - theta_bcd;

        cross_product(u1, y1, z);
        _normalize(z);

        cross_product(z, y1, f2);
        pos = _add(
                    _add(c,const_multiply(cos(rotate*3.142/180), y1)),
                    _add(const_multiply(sin(rotate*3.142/180), f2),
                        (const_multiply(1.0-cos(rotate*3.142/180),
                        const_multiply(dot_product(z, y1), z)))));

        _subtract(pos, c, c4);
        pos_BL = _add(const_multiply(bond_len/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), c);

        return pos_BL;
      }

V2 class2( V1 pos_CA, V1 pos_CG, V1 pos_CB, float bond_len)
      {
        V1 s, m, n, hb1, hb2, cdiv, ans, c3, c4, norm, norm2, point_hb1, point_hb2;

        m =  const_multiply(0.5,_add(pos_CA, pos_CG));
        _subtract(pos_CB, m, s);
        n = _add(m, const_multiply(2.0, s));
        _subtract(pos_CB, m, norm);
        _subtract(pos_CA, m, norm2);
        cross_product(norm2, norm, ans);

        cdiv = const_divide(sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]), ans);
        _subtract(n , cdiv, hb1);
        hb2 = _add(n, cdiv);

        if (bond_len !=0.0)
            {
              _subtract(hb1, pos_CB, c4);
              _subtract(hb2, pos_CB, c3);
              _add(const_multiply(bond_len/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), pos_CB);
                point_hb1 = _add(const_multiply(1.01/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), pos_CB);
                point_hb2 = _add(const_multiply(1.01/sqrt(c3[0]*c3[0]+c3[1]*c3[1]+c3[2]*c3[2]), c3), pos_CB);

                return {point_hb1, point_hb2};
            }

        _subtract(hb1, pos_CB, c4);
        _subtract(hb2, pos_CB, c3);
        _add(const_multiply(bond_len/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), pos_CB);
        point_hb1 = _add(const_multiply(1.09/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), pos_CB);
        point_hb2 = _add(const_multiply(1.09/sqrt(c3[0]*c3[0]+c3[1]*c3[1]+c3[2]*c3[2]), c3), pos_CB);

        return {point_hb1, point_hb2};
      }

V1 class3(V1 pos_CB, V1 pos_N, V1 pos_CA)
      {
        V1 m, n, s, ha,c4, ans, cdiv, norm, norm2, point_ha;

        m =  const_multiply(0.5,_add(pos_CB, pos_N));
        _subtract(pos_CA, m, s);
        n = _add(m, const_multiply(2.0, s));

        _subtract(pos_CA, m, norm);
        _subtract(pos_CB, m, norm2);
        cross_product(norm2, norm, ans);

        cdiv = const_divide(sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]), ans);
        ha   = _add(n, cdiv);

        _subtract(ha, pos_CA, c4);
        point_ha = _add(const_multiply(1.09/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), pos_CA);

        return point_ha;
      }

V1 class4(V1 pos_NE, V1 pos_CZ, V1 pos_NH1)
      {
        V1 p1, p2, p3, s1, s2, pos_hh11, point_1hh1;

        _subtract(pos_NE, pos_CZ, p1);
        _subtract(pos_NH1, pos_CZ, p2);
        p3 = _add(p1, p2);
        pos_hh11   = _add(p3, pos_NH1);

        _subtract(pos_hh11, pos_NH1, s1);
        _subtract(pos_hh11, pos_NH1, s2);

        point_1hh1 = _add(const_multiply(1.01/sqrt(s1[0]*s1[0]+s1[1]*s1[1]+s1[2]*s1[2]), s2), pos_NH1);

        return point_1hh1;
      }

V1 class5(V1 pos_CD, V1 pos_NE, V1 pos_CZ, float bond_len)
      {
        V1 p1, p2, p3, s1, s2, c4, pos, z, u1, u2, y1, f2, pos_he, pos_BL;

        _subtract(pos_NE, pos_CD, p1);
        _subtract(pos_NE, pos_CZ, p2);
        p3     = _add(p1, p2);
        pos_he = _add(p3, pos_NE);

        _subtract(pos_CZ, pos_NE, u1);
        _subtract(pos_he, pos_NE, y1);
        _subtract(pos_CD, pos_NE, u2);

        float mag_y1  = sqrt((y1[0]*y1[0])+(y1[1]*y1[1])+(y1[2]*y1[2]));
        float mag_u1  = sqrt((u1[0]*u1[0])+(u1[1]*u1[1])+(u1[2]*u1[2]));
        float mag_u2  = sqrt((u2[0]*u2[0])+(u2[1]*u2[1])+(u2[2]*u2[2]));

        float theta_bcd  = (acos(dot_product(u1,y1)/(mag_u1*mag_y1)))*180./3.142;
        float theta_bcd2 = acos(dot_product(u2,y1)/(mag_u2*mag_y1))*180/3.142;

        float theta     = (theta_bcd + theta_bcd2)/2;
        float rotate    = theta - theta_bcd;

        cross_product(u1, y1, z);
        _normalize(z);

        cross_product(z, y1, f2);
        pos = _add(
                      _add(pos_NE,const_multiply(cos(rotate*3.142/180.), y1)),
                      _add(const_multiply(sin(rotate*3.142/180.), f2),
                          (const_multiply(1.0-cos(rotate*3.142/180.),
                          const_multiply(dot_product(z, y1), z)))));

          _subtract(pos, pos_NE, c4);
        pos_BL = _add(const_multiply(bond_len/sqrt(c4[0]*c4[0]+c4[1]*c4[1]+c4[2]*c4[2]), c4), pos_NE);

        return pos_BL;
      }

V1 class6_Ser(V1 pos_OG, V1 pos_CB, V1 HB1)
  {
      V1 og_cb_vector  = {pos_CB[0]-pos_OG[0], pos_CB[1]-pos_OG[1], pos_CB[2]-pos_OG[2]};
      V1 cb_hb1_vector = {HB1[0]-pos_CB[0], HB1[1]-pos_CB[1], HB1[2]-pos_CB[2]};

      float rotation_amount = -240.2*(3.142/180.0);

      V1 rotated = rotMatrix(og_cb_vector, cb_hb1_vector, rotation_amount);
      V1 pos_h   = {rotated[0]+pos_OG[0], rotated[1]+pos_OG[1], rotated[2]+pos_OG[2]};

      V1 vec1    = {pos_h[0]-pos_OG[0], pos_h[1]-pos_OG[1], pos_h[2]-pos_OG[2]};

      float k    = 0.96/sqrt(((pos_h[0]-pos_OG[0])*(pos_h[0]-pos_OG[0])+
                      (pos_h[1]-pos_OG[1])*(pos_h[1]-pos_OG[1])+
                      (pos_h[2]-pos_OG[2])*(pos_h[2]-pos_OG[2])));

      V1 vec3    = {vec1[0]*k, vec1[1]*k, vec1[2]*k};

      V1 pos_hg  = {vec3[0]+pos_OG[0], vec3[1]+pos_OG[1], vec3[2]+pos_OG[2]};

      return pos_hg;

  }

  V1 class6_Thr(V1 pos_OG1, V1 pos_CB, V1 pos_CG2)
  {
      V1 og1_cb_vector  = {pos_CB[0]-pos_OG1[0], pos_CB[1]-pos_OG1[1], pos_CB[2]-pos_OG1[2]};
      V1 cb_cg2_vector = {pos_CG2[0]-pos_CB[0], pos_CG2[1]-pos_CB[1], pos_CG2[2]-pos_CB[2]};

      float rotation_amount = -243.2*(3.142/180.0);

      V1 rotated = rotMatrix(og1_cb_vector, cb_cg2_vector, rotation_amount);
      V1 pos_h   = {rotated[0]+pos_OG1[0], rotated[1]+pos_OG1[1], rotated[2]+pos_OG1[2]};

      V1 vec1    = {pos_h[0]-pos_OG1[0], pos_h[1]-pos_OG1[1], pos_h[2]-pos_OG1[2]};

      float k    = 0.96/sqrt(((pos_h[0]-pos_OG1[0])*(pos_h[0]-pos_OG1[0])+
                      (pos_h[1]-pos_OG1[1])*(pos_h[1]-pos_OG1[1])+
                      (pos_h[2]-pos_OG1[2])*(pos_h[2]-pos_OG1[2])));

      V1 vec3    = {vec1[0]*k, vec1[1]*k, vec1[2]*k};

      V1 pos_hg  = {vec3[0]+pos_OG1[0], vec3[1]+pos_OG1[1], vec3[2]+pos_OG1[2]};

      return pos_hg;

  }

  V1 class6_Tyr(V1 pos_OH, V1 pos_CZ, V1 pos_CE2)
  {
      V1 oh_cz_vector  = {pos_CZ[0]-pos_OH[0], pos_CZ[1]-pos_OH[1], pos_CZ[2]-pos_OH[2]};
      V1 cz_ce2_vector = {pos_CE2[0]-pos_CZ[0], pos_CE2[1]-pos_CZ[1], pos_CE2[2]-pos_CZ[2]};

      float rotation_amount = -220.2*(3.142/180.0);

      V1 rotated = rotMatrix(oh_cz_vector, cz_ce2_vector, rotation_amount);
      V1 pos_h   = {rotated[0]+pos_OH[0], rotated[1]+pos_OH[1], rotated[2]+pos_OH[2]};

      V1 vec1    = {pos_h[0]-pos_OH[0], pos_h[1]-pos_OH[1], pos_h[2]-pos_OH[2]};

      float k    = 0.96/sqrt(((pos_h[0]-pos_OH[0])*(pos_h[0]-pos_OH[0])+
                      (pos_h[1]-pos_OH[1])*(pos_h[1]-pos_OH[1])+
                      (pos_h[2]-pos_OH[2])*(pos_h[2]-pos_OH[2])));

      V1 vec3    = {vec1[0]*k, vec1[1]*k, vec1[2]*k};

      V1 pos_hg  = {vec3[0]+pos_OH[0], vec3[1]+pos_OH[1], vec3[2]+pos_OH[2]};

      return pos_hg;

  }

    V1 class6_Cys(V1 pos_SG, V1 pos_CB, V1 pos_CA)
  {
      V1 sg_cbvector  = {pos_CB[0]-pos_SG[0], pos_CB[1]-pos_SG[1], pos_CB[2]-pos_SG[2]};
      V1 cb_ca_vector = {pos_CA[0]-pos_CB[0], pos_CA[1]-pos_CB[1], pos_CA[2]-pos_CB[2]};

      float rotation_amount = -243.2*(3.142/180.0);

      V1 rotated = rotMatrix(sg_cbvector, cb_ca_vector, rotation_amount);
      V1 pos_h   = {rotated[0]+pos_SG[0], rotated[1]+pos_SG[1], rotated[2]+pos_SG[2]};

      V1 vec1    = {pos_h[0]-pos_SG[0], pos_h[1]-pos_SG[1], pos_h[2]-pos_SG[2]};

      float k    = 1.34/sqrt(((pos_h[0]-pos_SG[0])*(pos_h[0]-pos_SG[0])+
                      (pos_h[1]-pos_SG[1])*(pos_h[1]-pos_SG[1])+
                      (pos_h[2]-pos_SG[2])*(pos_h[2]-pos_SG[2])));

      V1 vec3    = {vec1[0]*k, vec1[1]*k, vec1[2]*k};

      V1 pos_hg  = {vec3[0]+pos_SG[0], vec3[1]+pos_SG[1], vec3[2]+pos_SG[2]};

      return pos_hg;

  }