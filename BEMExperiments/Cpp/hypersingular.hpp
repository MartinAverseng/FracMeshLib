#ifndef HYPERSINGULAR_HPP
#define HYPERSINGULAR_HPP

#include <cmath>
#include "quadrule.hpp"

double  PI      = 3.141592653589793;
double* rule[4] = {rule0+1,rule1+1,rule2+1,rule3+1};
int     nq[4]   = {(int)*rule0, (int)*rule1, (int)*rule2, (int)*rule3};

class R3{

private:

  double x[3];

public:

  R3(const double& x0=0.,const double& x1=0.,const double& x2=0.){
    x[0]=x0;x[1]=x1;x[2]=x2;}

  R3(const R3& p){x[0]=p.x[0];x[1]=p.x[1];x[2]=p.x[2];}

  R3& operator=(const R3& p){if(this!=&p){
      x[0]=p.x[0];x[1]=p.x[1];x[2]=p.x[2];}
    return *this;}

  double& operator[](const int& j){return x[j];}
  const double& operator[](const int& j) const {return x[j];}

  R3& operator+=(const R3& u){
    x[0]+=u[0];x[1]+=u[1];x[2]+=u[2]; return *this;}
  friend R3 operator+(R3 u, const R3& v){return u+=v;}

  R3& operator-=(const R3& u){
    x[0]-=u[0];x[1]-=u[1];x[2]-=u[2]; return *this;}
  friend R3 operator-(R3 u, const R3& v){return u-=v;}

  R3& operator*=(const double& a){
    x[0]*=a; x[1]*=a; x[2]*=a; return *this;}
  friend R3 operator*(const double& a, R3 u){return u*=a;}

  friend double operator,(const R3& u, const R3& v){
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];}

  friend double norm(const R3& u){return sqrt( (u,u) );}

  friend void normalize(R3& u){double a=norm(u); u*=1./a;}

};

R3 vprod(const R3& u, const R3& v){
  return R3(u[1]*v[2]-u[2]*v[1],
	    u[2]*v[0]-u[0]*v[2],
	    u[0]*v[1]-u[1]*v[0]);}


double  SingInt(double* vtx, int* num){

  int Iy[] = {num[0],num[1],num[2]};
  int Ix[] = {num[3],num[4],num[5]};

  int r  = 0;
  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      if(Ix[j]==Iy[k]){
	std::swap(Ix[j],Ix[r]);
	std::swap(Iy[k],Iy[r]);
	r++;
      }
    }
  }

  R3 x[3], y[3];
  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      x[j][k]=vtx[3*Ix[j]+k];
      y[j][k]=vtx[3*Iy[j]+k];
    }
  }

  double* qr = rule[r];
  double inter=0; R3 Xq,Yq;
  for(int q=0; q<nq[r]; q+=5){

    const double& u0 = qr[q  ];
    const double& u1 = qr[q+1];
    const double& v0 = qr[q+2];
    const double& v1 = qr[q+3];
    const double& wq = qr[q+4];

    Xq = x[0]+u0*(x[1]-x[0])+u1*(x[2]-x[0]);
    Yq = y[0]+v0*(y[1]-y[0])+v1*(y[2]-y[0]);

    inter += wq/norm(Xq-Yq);
  }

  return inter;

}

//#####################################################//
//   Assembly routine for Laplace
//   hypersingular integral operator
//
//   Input:
//
//   - vtx: A table of size 18=3x6 double
//          which contains the 3D coordinates
//          of 6 points x0, x1, ..., x5.
//          If (xj0, xj1, xj2) are the coordinates
//          of point xj, the the table vtx contains
//          vtx = [ x00, x01, x02, x10, ...., x52 ]
//
//   - tri: A table of size 6=2x3 int which are the
//          the indices j of the vertices of two triangles
//          T0=[v00, v01, v02] and T1=[v10, v11, v12]
//          where the vjk are one of the 6 points of vtx.
//          T0 is the source triangle and T1 is the target
//          triangle. The convention is that:
//
//            v0k = xj where j=tri[k]   for k in {0,1,2}
//            v1k = xj where j=tri[k+3] for k in {0,1,2}
//
//          The normal vector considered on triangle Tp, p=0,1
//          is given by Np = (vp1-vp0)x(vp2-vp1).
//
//   Output:
//
//   - res: A table of size 9 = 3x3 modelling rwo-wise the 3x3
//          elementary matrix A obtained through this assembly
//          routine. The correspondance is the following:
//
//          A = [ res[0], res[1], res[2],
//                res[3], res[4], res[5],
//                res[6], res[7], res[8] ]
//
//##############################################################//

void HsOp(double* vtx, int* tri, double* res){

  int* Jx = tri+3;
  int* Jy = tri;
  R3 x[3], y[3];
  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      x[j][k]=vtx[3*Jx[j]+k];
      y[j][k]=vtx[3*Jy[j]+k];
    }
  }

  R3 Nx,Ny;
  Nx = vprod(x[1]-x[0],x[2]-x[1]);
  Ny = vprod(y[1]-y[0],y[2]-y[1]);
  double dx = norm(Nx); normalize(Nx);
  double dy = norm(Ny); normalize(Ny);

  R3 nx[3], ny[3];
  for(int j=0; j<3; ++j){
    nx[(j+2)%3]=vprod(x[(j+1)%3]-x[j],Nx);
    ny[(j+2)%3]=vprod(y[(j+1)%3]-y[j],Ny);
  }

  R3 curlx[3], curly[3];
  for(int j=0; j<3; ++j){
    curlx[j] = (1./(x[(j+1)%3]-x[j],nx[j]))*vprod(nx[j],Nx);
    curly[j] = (1./(y[(j+1)%3]-y[j],ny[j]))*vprod(ny[j],Ny);
  }

  double inter = SingInt(vtx,tri)*dx*dy/(4.*PI);

  for(int j=0; j<3; ++j){
    for(int k=0; k<3; ++k){
      res[3*j+k] = inter*(curlx[j],curly[k]);
    }
  }

}

#endif
