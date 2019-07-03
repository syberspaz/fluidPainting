PVector K_QuinacridoneRose = new PVector(0.22, 1.47, 0.57);
PVector S_QuinacridoneRose = new PVector(0.05, 0.003, 0.03);
PVector K_FrenchUltramarine = new PVector(0.86, 0.86, 0.06);
PVector S_FrenchUltramarine = new PVector(0.005, 0.005, 0.09);
PVector K_CeruleanBlue = new PVector(1.52, 0.32, 0.25);
PVector S_CeruleanBlue = new PVector(0.06, 0.26, 0.40);
PVector K_HookersGreen = new PVector(1.62, 0.61, 1.64);
PVector S_HookersGreen = new PVector(0.01, 0.012, 0.003);
PVector K_HansaYellow = new PVector(0.06, 0.21, 1.78);
PVector S_HansaYellow = new PVector(0.50, 0.88, 0.009);
PVector K_CadmiumYellow = new PVector(0.10,0.36,3.45);
PVector S_CadmiumYellow = new PVector(0.97,0.65,0.007);
PVector K_BrilliantOrange = new PVector(0.13,0.81,3.45);
PVector S_BrilliantOrange = new PVector(0.005,0.009,0.007);

PVector K_CadmiumRed = new PVector(0.14, 1.08, 1.68);
PVector S_CadmiumRed = new PVector(0.77, 0.015, 0.018);
PVector K_IndianRed = new PVector(0.46, 1.07, 1.50);
PVector S_IndianRed = new PVector(1.28, 0.38, 0.21);
PVector K_InterferenceLilac = new PVector(0.08, 0.11, 0.07);
PVector S_InterferenceLilac = new PVector(1.25, 0.42, 1.43);
PVector K_PhthaloGreen = new PVector(1.55, 0.47, 0.63);
PVector S_PhthaloGreen = new PVector(0.01,0.05,0.035);
PVector K_BurntUmber = new PVector(0.74, 1.54, 2.10);
PVector S_BurntUmber = new PVector(0.09, 0.09, 0.004);

final int NUMPIGMENTS = 10;

class Fluid2 {
    int size;
    float dt;
    float diff;
    float visc;
    PVector[] pigmentsK;
    PVector[] pigmentsS;
    
    /* temp storage for fluid sim */
    float[] s0;
    float[] s1;
    float[] s2;
    
    /* density of the pigments used */
    float[] densityP0;
    float[] densityP1;  
    float[] densityP2;
    
    /* velocity */
    float[] Vx;
    float[] Vy;

    float[] Vx0;
    float[] Vy0;
    
    /* constructor */
    Fluid2(float dt, float diffusion, float viscosity) {
   
      this.size = N;
      this.dt = dt;
      this.diff = diffusion;
      this.visc = viscosity;
      
      this.pigmentsK = new PVector[10];
      this.pigmentsS = new PVector[10];
      this.s0 = new float[N*N];
      this.s1 = new float[N*N];
      this.s2 = new float[N*N];
  
      this.densityP0 = new float[N*N];
      this.densityP1 = new float[N*N];
      this.densityP2 = new float[N*N];
      
      this.Vx = new float[N*N];
      this.Vy = new float[N*N];
      
      this.Vx0 = new float[N*N];
      this.Vy0 = new float[N*N];
      this.setPigment(0,K_CadmiumYellow, S_CadmiumYellow);
      //this.setPigment(1,K_CeruleanBlue, S_CeruleanBlue);
      this.setPigment(2, K_BrilliantOrange, S_BrilliantOrange);
      this.setPigment(1, K_PhthaloGreen, S_PhthaloGreen);
    }
    void setPigment(int ind, PVector K, PVector S){
      if(ind < 0) return;
      if(ind >= NUMPIGMENTS) return;
      this.pigmentsK[ind] = new PVector();
      this.pigmentsS[ind] = new PVector();
      this.pigmentsK[ind].set(K);
      this.pigmentsS[ind].set(S);
    }
    
    PVector compositeKM(int i, int j){
       /* composite the layers of pigment */
       float c0 = this.densityP0[IX(i,j)];
       float c1 = this.densityP1[IX(i,j)];
       float c2 = this.densityP2[IX(i,j)];
       
       /* get R0, T0 for each layer of pigment */
       PVector[] RT0 = KM(this.pigmentsK[0],this.pigmentsS[0],c0);
       PVector[] RT1 = KM(this.pigmentsK[1],this.pigmentsS[1],c1);
       PVector[] RT2 = KM(this.pigmentsK[2],this.pigmentsS[2],c2);
       
       /* layer0 + layer 1 */   
       PVector[] L01 = KMComposite(RT0[0], RT0[1], RT1[0], RT1[1]);
       /* (layer0 + layer1) + layer 2 */
       PVector[] res = KMComposite(L01[0], L01[1], RT2[0], RT2[1]);
       PVector col = new PVector();
       col.set(res[0]);
       col.add(res[1]);
       return col; //<>//
    }

    void renderD(){
      PVector col = new PVector(100,100,100);
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
           float x = i*SCALE;
           float y = j*SCALE;
           
           col = compositeKM(i,j);
// temp for debugging right now
          // col.x = this.densityP0[IX(i,j)];
          // col.y = this.densityP1[IX(i,j)];
          // col.z = this.densityP2[IX(i,j)];
           
           //float r = this.densityR[IX(i,j)];
           //float Y = this.densityG[IX(i,j)];
           //float b = this.densityB[IX(i,j)];
           //PVector rgb = rybToRgb(r,Y,b);
                  
           fill(col.x*255,col.y*255,col.z*255);
           noStroke();
           square(x,y,SCALE);
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
    void addDensityKM(int x, int y, int pigmentInd, float concentration){
        if(pigmentInd ==0) this.densityP0[IX(x,y)] += concentration;
        if(pigmentInd ==1) this.densityP1[IX(x,y)] += concentration;
        if(pigmentInd ==2) this.densityP2[IX(x,y)] += concentration;

    }
    
    
    //void addDensity(int x, int y, PVector rgb){
    //  PVector ryb = rgbToRyb(rgb.x,rgb.y,rgb.z);

    //  this.densityR[IX(x,y)] += ryb.x;
    //  this.densityG[IX(x,y)] += ryb.y;
    //  this.densityB[IX(x,y)] += ryb.z;
      
    //}
    
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
      float[] s0       = this.s0;
      float[] s1       = this.s1;
      float[] s2       = this.s2;
       
      float[] densityP0 = this.densityP0;
      float[] densityP1 = this.densityP1;
      float[] densityP2 = this.densityP2;

/* vstep - velocity field */
      /* add external forces */
      
      /* transport */
      
      /* diffuse */
      
      /* project */


/* sstep - scalar field */

      /* add forces */
      /*transport */
      /* diffuse */ 
      /* dissipate */





      //poisson2D(s0,densityP0, 1,1);
      //copyArray(s0,densityP0);
 //     diffuse(1, Vx0, Vx, visc,dt);
  //    diffuse(2, Vy0, Vy, visc,dt);
      
   //   project(Vx0, Vy0, Vx, Vy);
      
  //    advect(1, Vx, Vx0, Vx0, Vy0, dt);
  //    advect(2, Vy, Vy0, Vx0, Vy0, dt);
      
  //    project(Vx, Vy,  Vx0, Vy0 );
//     /* transfer */
    //  transferPigment(densityP1, densityP0, 0.99);
   
      
  //    diffuse(0, s0, densityP0, diff, dt);
  //    diffuse(0, s1, densityP1, diff, dt);
      //diffuse(0, s2, densityP2, diff, dt);
      
 //     advect(0, densityP0, s0, Vx, Vy,dt);
 //     advect(0, densityP1, s1, Vx, Vy,dt);
    //  transferPigment(densityP2, densityP1, 0.98);
      //CahnHilliard(densityP0, s0, 1,1 , dt, 0.1, densityP1);
      //advect(0, densityP2, s2, Vx, Vy,dt*2); 
    
     
    }
}
