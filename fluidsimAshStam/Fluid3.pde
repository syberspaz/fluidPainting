
class FluidLayer{
  
  /* vector field for the layer */
  float[] Vx;
  float[] Vy;
  
  float[] Vx0;
  float[] Vy0;
  
  /* external forces field */
  float[] Fx;
  float[] Fy;
  
  /* temp */
  float[] s0;
  
  /* scalar field of the layer */
  float[] density;
  
  /* K,S values of the pigment */  
  PVector[] pigmentK;
  PVector[] pigmentS;
  PVector[] currentColor;
  
  float viscosity;
  float diffusion;
  
  boolean isReady;
  
  FluidLayer(float diff, float visc){
    isReady = false;
    Vx  = new float[N*N];
    Vy  = new float[N*N];
    Vx0 = new float[N*N];
    Vy0 = new float[N*N];
    Fx  = new float[N*N];
    Fy  = new float[N*N];
    s0  = new float[N*N];
    density = new float[N*N];
    currentColor = new PVector[N*N];
    for(int i=0;i<N*N;i++) currentColor[i] = new PVector(0,0,0);
   
    pigmentK = new PVector[N*N];
    pigmentS = new PVector[N*N];
    for(int i=0;i<N*N;i++){
      pigmentK[i] = new PVector(0.06, 0.21, 1.78);
      pigmentS[i] = new PVector(0.50, 0.88, 0.009);
    }
 
    viscosity = visc;
    diffusion = diff;
    init();
  }
  
  void init(){
    /* clear layer of all data */
    setField(density,0.01);
    setField(Vx,0);
    setField(Vy,0);
    setField(Vx0,0);
    setField(Vy0,0);
    setField(Fx,0);
    setField(Fy,0);
    setField(s0,0);
  }
  void clear(){
    /* clear layer of all data */
    setField(density,0);
    setField(Vx,0);
    setField(Vy,0);
    setField(Fx,0);
    setField(Fy,0);
  }
  
  void setViscosity(float v){ viscosity = v; }
  void setDiffusion(float d){ diffusion = d; }
  
  void setPigment(PVector K, PVector S){
   // this.pigmentK.set(K);
    //this.pigmentS.set(S);
  }
  
  void addVelocity(int x, int y, float amountX, float amountY) {
    Vx[IX(x,y)] += amountX;
    Vy[IX(x,y)] += amountY;
  }
    
  void addForce(int x, int y, float fx, float fy){
    Fx[IX(x,y)] += fx;
    Fy[IX(x,y)] += fy;
  }
  
  
  
  
  void addPaint(int x, int y, PVector K, PVector S, float conc){
    /* add paint to layer */
    int idx = IX(x,y);
    float d = this.density[idx];
    float totalD =  d + conc;
    
    this.density[idx] = totalD;
    float r0 = d/totalD;
    float r1 = conc/totalD;

    /* set pigments to proper volume ratios to maintain a constant volume per cell */
    pigmentK[idx].x = pigmentK[idx].x*r0 + r1*K.x;
    pigmentK[idx].y = pigmentK[idx].y*r0 + r1*K.y;
    pigmentK[idx].z = pigmentK[idx].z*r0 + r1*K.z;

    pigmentS[idx].x = pigmentS[idx].x*r0 + r1*S.x;
    pigmentS[idx].y = pigmentS[idx].y*r0 + r1*S.y;
    pigmentS[idx].z = pigmentS[idx].z*r0 + r1*S.z;

//    pigmentK[IX(x,y)].add(K);
 //   pigmentS[IX(x,y)].add(S);
    
   // float c0 = this.density[IX(x,y)];
   // float c1 = conc;
   // float total = c0 + c1;
   // float r0 = c0/total;
   // float r1 = c1/total;
   
   ///* composite color for drawing */
   // PVector[] RT0 = KM(mixedPigmentK, mixedPigmentS, r0);
   // PVector[] RT1 = KM(K,S,r1);
   // PVector[] comp = KMComposite(RT0[0], RT0[1], RT1[0], RT1[1]);
   // currentColor[IX(x,y)].set(comp[0]);
   // currentColor[IX(x,y)].add(comp[1]);
    
   // this.pigmentK.set(K);
   // this.pigmentS.set(S);
   // this.density[IX(x,y)] += conc;
  }
  PVector getRGB(int i, int j){
    /* convert K,S into RGB */
     int idx = IX(i,j);
    float d = this.density[idx];
    
    
    PVector K = pigmentK[idx];
    PVector S = pigmentS[idx];
    PVector[] rt = KM(K,S,d);
    PVector col = new PVector(0,0,0);
    col.set(rt[0]);
    col.add(rt[1]);
    return col;
  }
  
  PVector[] getRT(int i, int j){
    int idx = IX(i,j);
    float d = this.density[idx];
    
    
    PVector K = pigmentK[idx];
    PVector S = pigmentS[idx];
    PVector[] rt = KM(K,S,d);
    return rt;
  
  }
  void mixWithLayer(FluidLayer layer){

    /* mix incoming color into this layer */
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        int idx = IX(i,j);
        PVector cK = pigmentK[idx];
        PVector cS = pigmentS[idx];
        float cd = this.density[idx];
       // PVector[] RT0 = KM(cK,cS, cd);
        
        PVector K = layer.pigmentK[idx];
        PVector S = layer.pigmentS[idx];
        float d = layer.density[idx];
        
        float totalD = cd + d;
        float r0 = cd/totalD;
        float r1 = d/totalD;
        pigmentK[idx].x = pigmentK[idx].x*r0 + r1*K.x;
        pigmentK[idx].y = pigmentK[idx].y*r0 + r1*K.y;
        pigmentK[idx].z = pigmentK[idx].z*r0 + r1*K.z;
    
        pigmentS[idx].x = pigmentS[idx].x*r0 + r1*S.x;
        pigmentS[idx].y = pigmentS[idx].y*r0 + r1*S.y;
        pigmentS[idx].z = pigmentS[idx].z*r0 + r1*S.z;

        //PVector[] RT1 = KM(K,S,d);
        
        ///* composite together */
        //PVector[] RT = KM_Composite(RT0[0], RT0[1], RT1[0], RT1[1]);
        
      }
    }
    
    
    /* final add concentration of pigment */
    for(int i=0;i<N*N;i++){
      this.density[i] += layer.density[i];
    }
    
    
  }
  
  /* scalar field step */
  void stepScalarField(float dt){
      /* transport */
      diffuse(0, s0, density, diffusion, dt);
 
      /* advect */ 
      advect(0, density, s0, Vx, Vy,dt);
  }
  
  /* velocity field step  */
  void stepVelocityField(float dt){
    /* transport */
      diffuse(1, Vx0, Vx, viscosity,dt);
      diffuse(2, Vy0, Vy, viscosity,dt);
     
      /* advect  */
      advect(1, Vx, Vx0, Vx0, Vy0, dt);
      advect(2, Vy, Vy0, Vx0, Vy0, dt);
  
      /* project */
      project(Vx, Vy,  Vx0, Vy0 );
  }
  
  /* step the field forward in time */
  void step(float dt){
      addForce1D(Vx, Fx, dt);
      addForce1D(Vy, Fy, dt);
    
      /* update the velocities */
      stepVelocityField(dt);

      /* sstep - scalar field */
      stepScalarField(dt);
      
      /* RE compute colors based ????? */
  }
  
}



class Fluid3 {
    int size;
    float dt;

    PVector[] pigmentsK;
    PVector[] pigmentsS;
    
    /* let's try a 2-layer thing */
    /* we paint into one layer and then mix into another */
    FluidLayer mixedLayer;
    FluidLayer paintLayer;
    
    
    /* constructor */
    Fluid3(float dt, float diffusion, float viscosity) {
   
      this.size = N;
      this.dt = dt;

      this.pigmentsK = new PVector[10];
      this.pigmentsS = new PVector[10];



//PVector K_IndianRed = new PVector(0.46, 1.07, 1.50);
//PVector S_IndianRed = new PVector(1.28, 0.38, 0.21);
//PVector K_InterferenceLilac = new PVector(0.08, 0.11, 0.07);
//PVector S_InterferenceLilac = new PVector(1.25, 0.42, 1.43);





      this.setPigment(1,K_CadmiumYellow, S_CadmiumYellow);
      this.setPigment(2, K_HansaYellow, S_HansaYellow);
      
      this.setPigment(3, K_HookersGreen, S_HookersGreen);
      this.setPigment(4, K_PhthaloGreen, S_PhthaloGreen);
      
      this.setPigment(5, K_BrilliantOrange, S_BrilliantOrange);
      this.setPigment(6, K_CadmiumRed, S_CadmiumRed);

      this.setPigment(7, K_QuinacridoneRose, S_QuinacridoneRose);
      
      this.setPigment(8, K_FrenchUltramarine, S_FrenchUltramarine);
      this.setPigment(9, K_CeruleanBlue, S_CeruleanBlue);
      
      
      this.setPigment(0, K_BurntUmber, S_BurntUmber);
     
      mixedLayer = new FluidLayer(diffusion, viscosity*10);
      paintLayer = new FluidLayer(diffusion, viscosity);
      paintLayer.isReady = true;
      
    }
    void setPigment(int ind, PVector K, PVector S){
      if(ind < 0) return;
      if(ind >= NUMPIGMENTS) return;
      this.pigmentsK[ind] = new PVector();
      this.pigmentsS[ind] = new PVector();
      this.pigmentsK[ind].set(K);
      this.pigmentsS[ind].set(S);
    }
    void mixLayers(){
      mixedLayer.mixWithLayer(paintLayer);
      paintLayer.clear();
    }
    
   
    void render(){
      PVector col = new PVector(0,0,0);
      PVector col2 = new PVector(0,0,0);
      PVector[] RT0;
      PVector[] RT1;
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
           float x = i*SCALE;
           float y = j*SCALE;
           
           RT0 = mixedLayer.getRT(i,j);
           RT1 = paintLayer.getRT(i,j);
         
          /* composite together for viewing */
           PVector[] RT = KMComposite(RT0[0], RT0[1], RT1[0], RT1[1]);
           col.set(RT[0]);
           col.add(RT[1]);
        
           fill(col.x*255,col.y*255,col.z*255);
           noStroke();
           square(x,y,SCALE);
        }
      } 
    }

    void setViscosity(float v){
      paintLayer.setViscosity(v);
    }
    void setDiffusion(float d){
      paintLayer.setDiffusion(d);
    }
    
   
    void addPaint(int x, int y, int pigmentInd, float concentration){
      
      paintLayer.addPaint(x,y, pigmentsK[pigmentInd],pigmentsS[pigmentInd],concentration);
    }
   
    void addVelocity(int x, int y, float amountX, float amountY) {
      paintLayer.addVelocity(x,y,amountX,amountY);
      mixedLayer.addVelocity(x,y,amountX*0.3,amountY*0.3);
    }
    
    void addForce(int x, int y, float fx, float fy){
      paintLayer.addForce(x,y,fx,fy);
    }
  
    
    void step(){
      
      paintLayer.step(dt);
      mixedLayer.step(dt);
      
    }
}
