
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
    
    /* K, S, are the surface properties of the input paint */
    void addDensityKM(int x, int y, PVector K, PVector S){
      PVector R0 = new PVector();
      R0.x = this.densityR[IX(x,y)];
      R0.y = this.densityG[IX(x,y)];
      R0.z = this.densityB[IX(x,y)];
 
      /* use KM here to composite colors */
      /* what is xx ? brightness? */
      float xx = 0.5;
      PVector[] rt = KM(K, S, xx);
      PVector R = rt[0];  // reflectance
      PVector T = rt[1];  // transmittance
      
      /* now we have the KM model of the incoming paint, we need to composite it with what is already there */
     // PVector[] comp = KMComposite(R0, T0, R,T);
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
