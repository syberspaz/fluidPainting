import java.lang.Math;

// Table of pigments https://www.shadertoy.com/view/XscSDN
// from Computer-Generated Watercolor. Cassidy et al.
// K is absortion. S is scattering
//vec3 K_QuinacridoneRose = vec3(0.22, 1.47, 0.57);
//vec3 S_QuinacridoneRose = vec3(0.05, 0.003, 0.03);
//vec3 K_FrenchUltramarine = vec3(0.86, 0.86, 0.06);
//vec3 S_FrenchUltramarine = vec3(0.005, 0.005, 0.09);
//vec3 K_CeruleanBlue = vec3(1.52, 0.32, 0.25);
//vec3 S_CeruleanBlue = vec3(0.06, 0.26, 0.40);
//vec3 K_HookersGreen = vec3(1.62, 0.61, 1.64);
//vec3 S_HookersGreen = vec3(0.01, 0.012, 0.003);
//vec3 K_HansaYellow = vec3(0.06, 0.21, 1.78);
//vec3 S_HansaYellow = vec3(0.50, 0.88, 0.009);
//vec3 K_CadmiumRed = vec3(0.14, 1.08, 1.68);
//vec3 S_CadmiumRed = vec3(0.77, 0.015, 0.018);
//vec3 K_IndianRed = vec3(0.46, 1.07, 1.50);
//vec3 S_IndianRed = vec3(1.28, 0.38, 0.21);
//vec3 K_InterferenceLilac = vec3(0.08, 0.11, 0.07);
//vec3 S_InterferenceLilac = vec3(1.25, 0.42, 1.43);

 PVector tmp = new PVector();
 PVector R = new PVector();
 PVector T = new PVector();
// PVector[] ret = new PVector[2];
PVector[] KMComposite(PVector R0,PVector T0, PVector R1,PVector T1){
  //PVector tmp = new PVector();
  tmp.x = 1.0/(1.0 - R0.x*R1.x);
  tmp.y = 1.0/(1.0 - R0.y*R1.y);
  tmp.z = 1.0/(1.0 - R0.z*R1.z);
 // PVector R = new PVector();
  R.x = R0.x + T0.x*T0.x*R1.x*tmp.x;
  R.y = R0.y + T0.y*T0.y*R1.y*tmp.y;
  R.z = R0.z + T0.z*T0.z*R1.z*tmp.z;
  
//  PVector T = new PVector();
  T.x = T0.x*T1.x*tmp.x;
  T.y = T0.y*T1.y*tmp.y;
  T.z = T0.z*T1.z*tmp.z;
  
  PVector[] ret = new PVector[2];
  ret[0] = new PVector();
  ret[1] = new PVector();
  ret[0].set(R);
  ret[1].set(T);
  return ret;
  
//  vec3 tmp = vec3(1.0) / (vec3(1.0) - R0 * R1);
//    R = R0 + T0 * T0 * R1 * tmp;
//    T = T0 * T1 * tmp;
}
PVector a = new PVector();
PVector bb = new PVector();
PVector bSx = new PVector();
PVector sbSX = new PVector();
PVector c = new PVector();
PVector[] KM(PVector K, PVector S, float x){
// output is R,T 
//  PVector a = new PVector();
  a.set(K);
  a.add(S); // K + S;
  a.x /= S.x;
  a.y /= S.y;
  a.z /= S.z;
  
  //PVector b = new PVector();
  bb.x = sqrt(a.x*a.x - 1.0);
  bb.y = sqrt(a.y*a.y - 1.0);
  bb.z = sqrt(a.z*a.z - 1.0);
  
  //PVector bSx = new PVector();
  bSx.x = bb.x*(S.x*x);
  bSx.y = bb.y*(S.y*x);
  bSx.z = bb.z*(S.z*x);
  
 // PVector sbSX = new PVector();
  sbSX.x = (float) Math.sinh(bSx.x);
  sbSX.y = (float) Math.sinh(bSx.y);
  sbSX.z = (float) Math.sinh(bSx.z);
  
 // PVector c = new PVector();
  c.x = a.x*sbSX.x + bb.x*(float)Math.cosh(bSx.x);
  c.y = a.y*sbSX.y + bb.y*(float)Math.cosh(bSx.y);
  c.z = a.z*sbSX.z + bb.z*(float)Math.cosh(bSx.z);
  
  //PVector R = new PVector();
  
  R.set(sbSX);
  R.x /= c.x;
  R.y /= c.y;
  R.z /= c.z;
  
 // PVector T = new PVector();
  T.set(bb);
  T.x /= c.x;
  T.y /= c.y;
  T.z /= c.z;
  
  PVector[] ret = new PVector[2];
  ret[0] = new PVector();
  ret[1] = new PVector();
  ret[0].set(R);
  ret[1].set(T);
  return ret;
  
  
  //PVector a = (K + S) / S;
 // PVector b = sqrt(a * a - PVector(1.0, 1.0, 1.0));
  //PVector bSx = b * S * vec3(x, x, x);
  //PVector sinh_bSx = sinh(bSx);
//  PVector c = a * sinh_bSx + b * cosh(bSx);
      
  //R = sinh_bSx / c;
 // T = b / c;
}




PVector rybToRgb(float r, float y, float b){
  float iW = min(r,y,b);
  float R = r - iW;
  float Y = y - iW;
  float B = b - iW;
  
  float mY = max(R,Y,B);
  float G = min(Y,B);
  
  Y -= G;
  B -= G;
  
  if(B > 0 && G > 0){
    B *= 2.0f;
    G *= 2.0f;
  }
  
  R += Y;
  G += Y;
  
  float mG = max(R,G,B);
  if(mG > 0){
    float N = mY / mG;
    R *= N;
    G *= N;
    B *= N;
  }
  
  R += iW;
  G += iW;
  B += iW;
  
  return new PVector(R,G,B);
}



PVector rgbToRyb (float r,float g,float b) {
  float iW = min(r,g,b);
  float R = r - iW;
  float G = g - iW;
  float B = b - iW;
  
  float mG = max(R,G,B);
  float iY = min(R,G);
  R -=iY;
  G -= iY;
  if(B >0 && G > 0){
    B /=2.f;
    G /=2.f;
  }
  
  iY += G;
  B += G;
  
  float mY = max(R,iY,B);
  if(mY >0){
    float N = mG/mY;
    R *= N;
    iY *= N;
    B *= N;
  }
  R += iW;
  iY += iW;
  B += iW;
  return new PVector(R,iY,B);
}

void setField(float[] src, float value){
  for(int i=0;i<N*N;i++) src[i] = value;
}

void copyArray(float[] src, float[] dst){
  for(int i=0;i<N*N;i++) dst[i] = src[i];
}

void convection_linear1D(float[] uOut, float[] uIn, float c, float dt, float dx){
      for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                  int idx = IX(i,j);
                  uOut[idx] = uIn[idx] - c*(dt/dx)*(uIn[IX(i,j)] - uIn[IX(i-1,j)]);
                }
      }
}

void convection_nonlinear1D(float[] uOut, float[] uIn, float dt, float dx){
      for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                  int idx = IX(i,j);
                  uOut[idx] = uIn[idx] - uIn[idx]*(dt/dx)*(uIn[IX(i,j)] - uIn[IX(i-1,j)]);
                }
      }
}


void diffusion_linear1D(float[] uOut, float[] uIn, float visc, float dt, float dx){
      for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                  int idx = IX(i,j);
                  uOut[idx] = uIn[idx] - (visc*dt/(dx*dx))*(uIn[IX(i+1,j)] - 2.0*uIn[IX(i,j)]+ uIn[IX(i-1,j)]);
                }
      }
}


void burger_linear1D(float[] uOut, float[] uIn, float visc, float dt, float dx){
      for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                  int idx = IX(i,j);
                  uOut[idx] = uIn[idx] - uIn[IX(i,j)]*(dt/dx)*(uIn[IX(i,j)] - uIn[IX(i-1,j)])+ (visc*dt/(dx*dx))*(uIn[IX(i+1,j)] - 2.0*uIn[IX(i,j)]+ uIn[IX(i-1,j)]);
                }
      }
}

void laplacian2D(float[] uOut, float[] uIn,  float dx, float dy){
      for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                  int idx = IX(i,j);
              
                  uOut[idx] = constrain((laplacian2d(uIn, i,j,dx,dy)),0,255);
                }
      }
}

void poisson2D(float[] uOut, float[] uIn,  float dx, float dy){
      for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                  int idx = IX(i,j);
                  float dy2 = dy*dy;
                  float dx2 = dx*dx;
                  float denom = 2*(dx2+dy2);
                  float bij = 0;//uIn[IX(i,j)];
                  uOut[idx] = ((uIn[IX(i+1,j)] + uIn[IX(i-1,j)])*dy2 + (uIn[IX(i,j+1)] + uIn[IX(i,j-1)])*dx2 - bij*dx2*dy2)/denom;
                }
      }
}


////////////////////////////
// FLUID SIMULATION UTILS
void set_bnd(int b, float[] x)
{
      for(int i = 1; i < N - 1; i++) {
          x[IX(i, 0 )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1  )];
          x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
      }
  
      for(int j = 1; j < N - 1; j++) {
          x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
          x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
      }
    x[IX(0, 0)]     = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)]   = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);
    x[IX(N-1, 0)]   = 0.5f * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);    
}

void lin_solve(int b, float[] x, float[] x0, float a, float c)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a*(    x[IX(i+1, j  )]
                                    +x[IX(i-1, j )]
                                    +x[IX(i  , j+1 )]
                                    +x[IX(i  , j-1 )]
                             
                           )) * cRecip;
                }
            }
        
        set_bnd(b, x);
    }
}


void diffuse (int b, float[] x, float[] x0, float diff, float dt)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

void computeForcesFromField(float[] s, float[] Fx, float[] Fy){
   int dx,dy;
   dx = dy = 3;
  for (int j = 2; j < N - 2; j++) {
            for (int i = 2; i < N - 2; i++) {
              /* gradient of scalar field */
              float gx, gy;
              //gx = (s[IX(i+1,j)] - s[IX(i-1,j)])/(2*dx);
             // s[IX(i+1,j)] - s[IX(i-1,j)])/(2*dx);
              
              for (int jj = j-dy; jj <= j+dy; jj++) {
                for (int ii = i-dx; ii <= i+dx; ii++) {
                  gx = (s[IX(ii+1,jj)] - s[IX(ii-1,jj)])/2;
                  gy = (s[IX(ii,jj+1)] - s[IX(ii,jj-1)])/2;
                  Fx[IX(ii,jj)] += gx*s[IX(ii,jj)];
                  Fy[IX(ii,jj)] += gy*s[IX(ii,jj)];
                }
              }
            }
  }
}

void project(float []velocX, float[] velocY, float[] p, float[] div)
{
 
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                         velocX[IX(i+1, j  )]
                        -velocX[IX(i-1, j   )]
                        +velocY[IX(i  , j+1 )]
                        -velocY[IX(i  , j-1  )]
                    )/N;
                p[IX(i, j)] = 0;
            }
        }
    
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);
    
  
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]-p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]-p[IX(i, j-1)]) * N;
              
            }
        }
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}
void transferPigment(float[] fromD, float[] toD, float rate){

  for(int i=0;i<N*N;i++){
    toD[i] += fromD[i]*(1-rate);  
    fromD[i] *= rate;
  }
}

/* 1D addforce/source to existing field */
void addForce1D(float[] v, float[] f, float dt){
  for(int i=0;i<N*N;i++){
    v[i] += f[i]*dt;
  }
}


float laplacian2d(float[] d, int x, int y, float dx, float dy){
  float current = d[IX(x,y)];
  float left    = d[IX(x-1,y)];
  float right   = d[IX(x+1,y)];
  float top     = d[IX(x,y+1)];
  float bottom  = d[IX(x,y-1)];

  float ret = (left + right - 2.0*current/(dx*dx)) + (top + bottom - 2.0*current)/(dy*dy);
  
  return ret;
}

void CahnHilliard_p1(float[] C, float[] intermediate, int x, int y, float dx, float dy, float gamma){
  float laplace = laplacian2d(C,x,y,dx,dy);
  
  float current = C[IX(x,y)];
  float current3 = current*current*current;
  float inter = current3 - current - gamma*laplace;
  
  intermediate[IX(x,y)] = inter; 
}

void CahnHilliard_p2(float[] C, float[] intermediate, int x, int y, float dx, float dy, float dt, float D){
  
  // phase 2
  float dC = D*laplacian2d(intermediate, x,y, dx,dy);
  float current = C[IX(x,y)];
  C[IX(x,y)] = current + dC*dt;
 
}

/* modifies incoming C density */
void CahnHilliard(float[] C, float[] intermediate, float dx, float dy, float dt, float gamma, float[] D){

    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
    for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
        for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            CahnHilliard_p1(C, intermediate, i,j, dx, dy, gamma);
            CahnHilliard_p2(C, intermediate, i,j, dx, dy,dt, D[IX(i,j)]);
        }
    }
}

void advect(int b, float[] d, float[] d0,  float[] velocX, float[] velocY,float dt)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floor(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floor(y);
                j1 = j0 + 1.0f; 
              
               
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
               
                
                int i0i = int(i0);
                int i1i = int(i1);
                int j0i = int(j0);
                int j1i = int(j1);
              

                d[IX(i, j)] = 
                    s0 * (t0*d0[IX(i0i, j0i)]+ t1*d0[IX(i0i, j1i)]) + 
                    s1 * (t0*d0[IX(i1i, j0i)]+ t1*d0[IX(i1i, j1i)]);
            }
        }
    
    set_bnd(b, d);
}
