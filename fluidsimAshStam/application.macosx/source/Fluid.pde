
float tt=0;
class Fluid {
    int size;
    float dt;
    float diff;
    float visc;
    
    float[] sR;
    float[] sG;
    float[] sB;
    float[] densityR;
    float[] densityG;  // testing: this will actually be Y, using RYB now
    float[] densityB;
    
    float[] Vx;
    float[] Vy;

    float[] Vx0;
    float[] Vy0;

void renderD(){
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
       float x = i*SCALE;
       float y = j*SCALE;
       
       float r = this.densityR[IX(i,j)];
       float Y = this.densityG[IX(i,j)];
       float b = this.densityB[IX(i,j)];
       PVector rgb = rybToRgb(r,Y,b);
              
       fill(rgb.x,rgb.y,rgb.z);
       noStroke();
       square(x,y,SCALE);
    }
  } 
}
void fadeD(){
 // for(int i=0;i<this.densityR.length;i++){
//    float d = density[i];
//    density[i] = constrain(d-0.01,0,255);
 // }
  
}


    Fluid(float dt, float diffusion, float viscosity) {
   
    this.size = N;
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;
    
    this.sR = new float[N*N];
    this.sG = new float[N*N];
    this.sB = new float[N*N];

    this.densityR = new float[N*N];
    this.densityG = new float[N*N];
    this.densityB = new float[N*N];
    
    this.Vx = new float[N*N];
    this.Vy = new float[N*N];
    
    this.Vx0 = new float[N*N];
    this.Vy0 = new float[N*N];
    
    }
    // add some weird random velocity thing
    void addVelocityField(){
      float cx,cy;
      cx = 0.5*width;
      cy = 0.5*height;
      
      for(int j=0;j<N;j++){
      for(int i=0;i<N;i++)
      {
        float vx = this.Vx[IX(i,j)];
        float vy = this.Vy[IX(i,j)];        
        float amtX = vx*(cx-i);
        float amtY = vy*(cy-j);
            this.Vx[IX(i,j)] += amtX;
            this.Vy[IX(i,j)] += amtY;

      }
    }
  }
  void setImage(PImage img){
      PVector rgb = new PVector();
 
    for(int j=0;j<N;j++){
    for(int i=0;i<N;i++){
      color c = img.get(i*SCALE,j*SCALE);
      rgb.set(red(c), green(c), blue(c));
      this.addDensity(i,j,rgb);
    }
    }
    
  }
  
    void setViscosity(float v){
    this.visc = v;
    }
    void setDiffusion(float d){
    this.diff = d;
    }
    
    void addDensity(int x, int y, PVector rgb){
      PVector ryb = rgbToRyb(rgb.x,rgb.y,rgb.z);
      
      this.densityR[IX(x,y)] += ryb.x;
      this.densityG[IX(x,y)] += ryb.y;
      this.densityB[IX(x,y)] += ryb.z;
      
    }
    
    void addVelocity(int x, int y, float amountX, float amountY) {
      int index = IX(x,y);
      this.Vx[index] += amountX;
      this.Vy[index] += amountY;
    }
    
    
    void step(){
      float visc     = this.visc;
      float diff     = this.diff;
      float dt       = this.dt;
      float[] Vx      = this.Vx;
      float[] Vy      = this.Vy;
      float[] Vx0     = this.Vx0;
      float[] Vy0     = this.Vy0;
      float[] sR       = this.sR;
      float[] sG       = this.sG;
      float[] sB       = this.sB;
       
      float[] densityR = this.densityR;
      float[] densityG = this.densityG;
      float[] densityB = this.densityB;

      
      diffuse(1, Vx0, Vx, visc,dt);
      diffuse(2, Vy0, Vy, visc,dt);
      
      project(Vx0, Vy0, Vx, Vy);
      
      advect(1, Vx, Vx0, Vx0, Vy0, dt);
      advect(2, Vy, Vy0, Vx0, Vy0, dt);
      
      project(Vx, Vy,  Vx0, Vy0 );
      
      diffuse(0, sR, densityR, diff, dt);
      diffuse(0, sG, densityG, diff, dt);
      diffuse(0, sB, densityB, diff, dt);
      
      advect(0, densityR, sR, Vx, Vy,dt);
      advect(0, densityG, sG, Vx, Vy,dt);
      advect(0, densityB, sB, Vx, Vy,dt);
     
    }
}
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
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N;
              
            }
        }
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
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
