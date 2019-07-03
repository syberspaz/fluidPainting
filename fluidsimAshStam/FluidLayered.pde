

class FluidLayered {
    float dt, viscosity, diffusion;
    /* fluid layers for paint */
    PaintPixel[] dry;
    PaintPixel[] wet;
    PaintPixel[] fluid;
    
    /* fluid layers for velocity */
    float[] Dry_Vx;
    float[] Dry_Vy;

    float[] Wet_Vx;
    float[] Wet_Vy;

    float[] Fluid_Vx;
    float[] Fluid_Vy;
    float[] Fluid_Vx0;
    float[] Fluid_Vy0;

    
    PaintPixel[] newLayer(int Nx, int Ny){
      
      PaintPixel[] layer = new PaintPixel[Nx*Ny];
      for(int i=0;i<Nx*Ny;i++){
        layer[i] = new PaintPixel();
      }
      return layer;
    
    }
    
    /* constructor */
    FluidLayered(float dt, float diffusion, float viscosity) {
       dry = newLayer(N,N);
       wet = newLayer(N,N);
       fluid = newLayer(N,N);
       
       Dry_Vx = new float[N*N];
       Dry_Vy = new float[N*N];
       Wet_Vx = new float[N*N];
       Wet_Vy = new float[N*N];
       Fluid_Vx = new float[N*N];
       Fluid_Vy = new float[N*N];
       Fluid_Vx0 = new float[N*N];
       Fluid_Vy0 = new float[N*N];
       
       this.dt = dt;
       this.viscosity = viscosity;
       this.diffusion = diffusion;
    }
   
    
    
    void render(){
      PVector col = new PVector(100,100,100);
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
           float x = i*SCALE;
           float y = j*SCALE;
           /* compute colors of each paint layer */
           /* mix all dry, wet, fluid colors together */
           PaintPixel p = fluid[IX(i,j)];
           col.x = col.y = col.z = p.volume;
          // col.x = Fluid_Vx[IX(i,j)]*255;//p.volume;
          // col.y = col.z = Fluid_Vy[IX(i,j)]*255;//col.x;
           fill(col.x*255,col.y*255,col.z*255);
           noStroke();
           square(x,y,SCALE);
        }
      } 
    }

    void setViscosity(float v){
      this.viscosity = v;
    }
    void setDiffusion(float d){
      this.diffusion = d;
    }
    
    
    void FLUID_DYNAMICS(float[] Vx, float[] Vy, float[] Vx0, float[] Vy0, float viscosity, float dt){
      diffuse(1, Vx0,Vx, viscosity, dt);
      diffuse(2, Vy0,Vy, viscosity, dt);
      project(Vx0, Vy0, Vx, Vy);
      advect(1,Vx, Vx0, Vx0, Vy0, dt);
      advect(2,Vy, Vy0, Vx0, Vy0, dt);
      project(Vx, Vy, Vx0, Vy0);
    }
    
    void AddPaint(int x, int y, Pigment p){
      fluid[IX(x,y)].addPaint(p);
    }
     void addVelocity(int x, int y, float amountX, float amountY) {
      int index = IX(x,y);
      this.Fluid_Vx[index] += amountX;
      this.Fluid_Vy[index] += amountY;
    }
    
    void AdvectPaint(PaintPixel[] layer, float[] Vx, float[] Vy, float dt){
      for(int i=0;i<N;i++){
          for(int j=0;j<N;j++){
              int idx = IX(i,j);
              PaintPixel p = layer[idx];
              float vx,vy;
              vx = Vx[idx];
              vy = Vy[idx];
              
              p.advectOut(vx,vy,dt);
              
              int iN = (j+1)<N? IX(i,j+1):0;
              int iS = (j-1)>=0? IX(i,j-1):0;
              int iE = (i-1)>=0? IX(i-1,j):0;
              int iW = (i+1)<N? IX(i+1,j):0;
              
              PaintPixel pX = (vx < 0)? layer[iE]: layer[iW];
              PaintPixel pY = (vy < 0)? layer[iS]: layer[iN];
                         
              p.advectIn(Math.abs(vx), Math.abs(vy), dt, pX, pY);
            
          }
      }
    
    
    }
    
    
    void step(){
      
      /* ok, so the algorithm is */
      /* Fluid Layer Dynamics */
      /* transfer paint from fluid layer to wet layer */
      
      /* Wet Layer Dynamics */
      /* transfer paint to dry layer */
      /* Dry Layer Dynamics (color mixing?, absorption, diffusion, etc) */
      
      
      /* fluid layer dynamics */
      FLUID_DYNAMICS(Fluid_Vx, Fluid_Vy, Fluid_Vx0, Fluid_Vy0, viscosity, dt);
      
      /* have to diffuse/advect the actual pigments now */
      AdvectPaint(fluid, Fluid_Vx, Fluid_Vy, dt);
      
      /* transfer pigment from fluid layer to wet layer */
      
      /* wet layer dynamics */      
     // FLUID_DYNAMICS(Wet_Vx, Wet_Vy, Fluid_Vx0, Fluid_Vy0, viscosity, dt);

      /* transfer to dry layer */
      /* dry layer dynamics */
      
      
      
    //  float visc     = this.visc;
    //  float diff     = this.diff;
    //  float dt       = this.dt;
    //  float[] Vx      = this.Vx;
    //  float[] Vy      = this.Vy;
    //  float[] Vx0     = this.Vx0;
    //  float[] Vy0     = this.Vy0;
    //  float[] s0       = this.s0;
    //  float[] s1       = this.s1;
    //  float[] s2       = this.s2;
       
    //  float[] densityP0 = this.densityP0;
    //  float[] densityP1 = this.densityP1;
    //  float[] densityP2 = this.densityP2;

      
    //  diffuse(1, Vx0, Vx, visc,dt);
    //  diffuse(2, Vy0, Vy, visc,dt);
      
    //  project(Vx0, Vy0, Vx, Vy);
      
    //  advect(1, Vx, Vx0, Vx0, Vy0, dt);
    //  advect(2, Vy, Vy0, Vx0, Vy0, dt);
      
    //  project(Vx, Vy,  Vx0, Vy0 );
    // /* transfer */
    ////  transferPigment(densityP1, densityP0, 0.99);
      
    //  diffuse(0, s0, densityP0, diff, dt);
    //  diffuse(0, s1, densityP1, diff, dt);
    //  diffuse(0, s2, densityP2, diff, dt);
      
    //  advect(0, densityP0, s0, Vx, Vy,dt);
    //  advect(0, densityP1, s1, Vx, Vy,dt);
    ////  transferPigment(densityP2, densityP1, 0.98);

    //  advect(0, densityP2, s2, Vx, Vy,dt*2);
     
     
     
    }
}
